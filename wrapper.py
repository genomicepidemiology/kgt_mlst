import os
import sys

path = '/home/people/malhal/cecilie/'
file_list = os.listdir(path)
file_list.sort()

for i in range(0, len(file_list), 2):
    name = file_list[i].split('/')[-1].split('.')[0]
    cmd = '/home/people/malhal/kgt_mlst/bin/kgt_mlst -i {} {} -db_dir {} -md 5 -o {}'\
        .format(path + '/' + file_list[i], path + '/' + file_list[i+1], '/home/people/malhal/kgt_mlst/kgt_mlst_db/', name)
    os.system(cmd)

