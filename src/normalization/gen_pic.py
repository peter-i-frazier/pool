import csv
import subprocess
def construct_stripes(img_folder, tmp_folder, dest_folder, lib_folder, img_list, ts, num_row, num_col):
    for img in img_list:
        subprocess.call(['convert', "{0}/{1}.png".format(img_folder, img), '-crop',
                         "1x{0}@".format(num_row), "{0}/row_{1}_%d.png".format(tmp_folder, img)])
        for row in range(num_row):
            subprocess.call(['convert', "{0}/row_{1}_{2}.png".format(tmp_folder, img, row),
                             '-crop', "{0}x1@".format(num_col),
                             "{0}/{1}_{2}_%d.png".format(tmp_folder, img, row)])
    for row in range(num_row):
        for col in range(num_col):
            target = "{0}/TS{1}_row_{2}_col_{3}.png".format(dest_folder, ts, row, col)
            subprocess.call(['convert', "{0}/ts_{1}.png".format(lib_folder, ts),
                             "{0}/r_{1}.png".format(lib_folder, row),
                             "{0}/c_{1}.png".format(lib_folder, col), '+append', target])
            for img in img_list:
                to_append = "{0}/{1}_{2}_{3}.png".format(tmp_folder, img, row, col)
                subprocess.call(['convert', target, to_append, '+append', target])

def rank_stripes(dest_folder, lib_folder, rank_table_folder, rank_table_name):
    f = open("{0}/{1}.csv".format(rank_table_folder, rank_table_name), 'rb')
    csv_f = csv.reader(f)
    for l, line in enumerate(csv_f):
        if l == 1:
            to_append = "{0}/TS{1}_row_{2}_col_{3}.png".format(lib_folder, line[1], line[3], line[4])
            target = "{0}/{1}.png".format(dest_folder, rank_table_name)
            subprocess.call(['convert', to_append, '-append', target])
        elif l > 1:
            to_append = "{0}/TS{1}_row_{2}_col_{3}.png".format(lib_folder, line[1], line[3], line[4])
            subprocess.call(['convert', target, to_append, '-append', target])
    f.close()
#name = "TS1_sfp_2/TS1_sfp_2"
#target = "TS1_sfp_2_rank.png"
#f = open('rank_1_2.csv', 'rb')
#csv_f = csv.reader(f)
#for l, line in enumerate(csv_f):
#    if (l > 0):
#        row = lookup_table[line[1][0]]
#        col = int(line[1][1:])
#        filename = name + "_" + str(row) + "_" + str(col-1) + ".png"
#        if (l == 1):
#            subprocess.call(['./init_append.sh', target, filename])
#        else:
#            subprocess.call(['./append.sh', target, filename])

#subprocess.call(['rm', '-rf', 'temp'])
#subprocess.call(['mkdir', 'temp'])
#construct_stripes('TS1', TS1_file, 9, 30)
#construct_stripes('TS2', TS2_file, 20, 30)
#construct_stripes('TS3', TS3_file, 19, 30)
#construct_stripes('TS4', TS4_file, 20, 30)
#construct_stripes('TS5', TS5_file, 20, 30)
#subprocess.call(['rm', '-rf', 'result'])
#subprocess.call(['mkdir', 'result'])
#ts_list = [1, 3, 4, 5]
#measure_list = [1, 2, 3, 4]
#for ts in ts_list:
#    for measure in measure_list:
#        append_stripes('result', ts, measure)
#append_stripes('result', 2, 1)
if __name__ == "__main__":
# construct stripes
    #subprocess.call(['rm', '-rf', 'tmp_img'])
    #subprocess.call(['mkdir', 'tmp_img'])
    #TS1 = ['TS1_sfp_1', 'TS1_sfp_2', 'TS1_AcpS_1', 'TS1_AcpS_2']
    #TS2 = ['TS2_sfp_1', 'TS2_sfp_2']
    #TS3 = ['TS3_sfp_1', 'TS3_sfp_2', 'TS3_AcpS_1', 'TS3_AcpS_2']
    #TS4 = ['TS4_sfp_1', 'TS4_sfp_2', 'TS4_AcpS_1', 'TS4_AcpS_2']
    #TS5 = ['TS5_sfp_1', 'TS5_sfp_2', 'TS5_AcpS_1', 'TS5_AcpS_2']
    #construct_stripes("org_img", "tmp_img", "lib_img", "lib_img", TS1, 1, 10, 30)
    #construct_stripes("org_img", "tmp_img", "lib_img", "lib_img", TS2, 2, 20, 30)
    #construct_stripes("org_img", "tmp_img", "lib_img", "lib_img", TS3, 3, 19, 30)
    #construct_stripes("org_img", "tmp_img", "lib_img", "lib_img", TS4, 4, 20, 30)
    #construct_stripes("org_img", "tmp_img", "lib_img", "lib_img", TS5, 5, 20, 30)
# rank stripes
    subprocess.call(['rm', '-rf', 'result_img'])
    subprocess.call(['mkdir', 'result_img'])
    rank_stripes("result_img", "lib_img", "rank_table", "rank_by_measure_1")
    rank_stripes("result_img", "lib_img", "rank_table", "rank_by_measure_2")
    rank_stripes("result_img", "lib_img", "rank_table", "rank_by_measure_3")
    rank_stripes("result_img", "lib_img", "rank_table", "rank_by_measure_4")
