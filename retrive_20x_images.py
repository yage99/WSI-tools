############################################################
# 
# Crop whole slide histopathological images(WSI, svs) to patch.
# All patches will have 20x magnification. WSIs with 40x magnification are down-
# sampled to 20x, and WSIs with 5x magnification are discard.
# Patches that are blank or without enough pathological tissue are filtered.
# 
# created by zhangya<zhangya998@gmail.com>
# 2021-11-07
#
############################################################
#!/usr/bin/python3
import os
import numpy as np
from TCGAMaxim import meta
import openslide
import tqdm
from multiprocessing import Pool, Array, Lock
import PIL
# from skimage import color
from TCGAMaxim import rgb2lab
import sys


tile_size = 1024
MAX_THREADS = os.cpu_count()
MIN_DENSITY = 0.35


pbar = None
def lab_density(image):
    pixel = np.array(image) / 255.0

    lab = rgb2lab(pixel[:, :, :3])

    pixel_count = lab.shape[0] * lab.shape[1]
    # the a channel should bigger than 10 (red)
    color_count = np.sum(np.logical_and(lab[:, :, 1] > 10,
                                        lab[:, :, 2] > -30))

    return color_count / float(pixel_count)


def calc_density(image, threshold=150):
    pixel = np.array(image)

    pixel_count = image.size[0] * image.size[1]
    color_count = np.sum(np.logical_or(
        pixel[:, :, 0] < threshold,
        pixel[:, :, 1] < threshold,
        pixel[:, :, 2] < threshold,
        ))
    transparent = np.sum(pixel[:, :, 3] < threshold)
    black = np.sum(np.sum(pixel[:, :, :3], axis=2) < 50)

    color_count = color_count - transparent - black
    # for i in range(image.size[0]):
    #     for j in range(image.size[1]):
    #         # transparent pixel don't count
    #         if pixel[i, j][3] < 125:
    #             continue
    #         # black pixel don't count
    #         if pixel[i, j][0] + pixel[i, j][1] + pixel[i, j][2] < 50:
    #             continue

    #         # red pixel has a minium
    #         if pixel[i, j][0] < threshold and pixel[i, j][0] > 50:
    #             color_count = color_count + 1
    #         # green pixel
    #         if pixel[i, j][1] < threshold:
    #             color_count = color_count + 1
    #         if pixel[i, j][2] < threshold:
    #             color_count = color_count + 1

    return color_count / pixel_count


def check_density(image):
    return lab_density(image) > MIN_DENSITY

def process_file(file, dest, dest_mag=20, pbar=None):
    global tqdm_slot, _lock

    slide = openslide.OpenSlide(file)

    mag = int(slide.properties['aperio.AppMag'])
    if mag < 15:
        return 'mag is %d, less than 15' % mag
    # downsample = slide.level_downsamples[lowest_level]
    downsample = mag / dest_mag + 1  # normalize to 5x
    best_level = slide.get_best_level_for_downsample(downsample)
    dimension = slide.level_dimensions[best_level]

    mag_dimension = slide.level_dimensions[0]
    level_mag = dimension[0] / mag_dimension[0] * mag
    level_tile = int(level_mag * tile_size / dest_mag)

    x = dimension[0] // level_tile
    y = dimension[1] // level_tile

    position = 0
    for i in range(len(tqdm_slot)):
        if tqdm_slot[i] == 1:
            position = i
            tqdm_slot[i] = 0
            break
    subbar = tqdm.tqdm(total=x*y, position=position + 1)
    for i in range(x):
        for j in range(y):
            subbar.update()
            img = slide.read_region(
                    (i * level_tile, j * level_tile),
                    best_level,
                    (level_tile, level_tile))
            img = img.resize((tile_size, tile_size), PIL.Image.NEAREST)
            if check_density(img):
                img.save("%s_%03d%03d.png" % (dest, i, j))

    tqdm_slot[position] = 1

    return best_level, mag, level_tile


def file_callback(*a):
    # pbar.write(str(a[0]))
    pbar.update()

tqdm_slot = Array('i', range(MAX_THREADS), lock=False)
_lock = Lock()
def init(l, t):
    global _lock, tqdm_slot
    _lock = l
    tqdm_slot = t

def process(meta_file, src_dir, dest_dir):
    global pbar, tqdm_slot

    meta_info = meta(meta_file)
    tasks_num = len(meta_info)
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    print('Total file number: %d' % tasks_num)
    pool = Pool(MAX_THREADS, initializer=init, initargs=(_lock, tqdm_slot,))
    pbar = tqdm.tqdm(total=tasks_num)

    for i in range(len(tqdm_slot)):
        tqdm_slot[i] = 1

    for id, file in meta_info:
        file_name = os.path.basename(file)
        obj = pool.apply_async(
                process_file,
                args=(os.path.join(src_dir, file),
                      os.path.join(dest_dir, file_name[:-4])
                      ),
                callback=file_callback
                )
        # process_file(os.path.join(src_dir, file),
        #     os.path.join(dest_dir, file_name[:-4]), tqdm_slot)

    pool.close()
    pool.join()
    pbar.close()

if __name__ == "__main__":
    if (len(sys.argv) > 1):
        process(sys.argv[1], 
            sys.argv[2],
            sys.argv[3]
            )
    else:
        print("""
Crop whole slide histopathological images(WSI, svs) to patch.
All patches will have 20x magnification. WSIs with 40x magnification are downsa-
mpled to 20x, and WSIs with 5x magnification are discard.
Patches that are blank or without enough pathological tissue are filtered.

For any question please contact <zhangya998@gmail.com>

Usage:
    retrive_20x_images <gdc_meta_txt_file> <gdc_download_folder> <output_folder>

    1. gdc_meta_txt_file
       get from TCGA gdc data portal, the same file for the TCGA official downl-
       oader gdc-client.
    2. gdc_download_folder
       the folder where gdc-client download all files
    3. output_folder
       The destination for selected patches.
        """)
