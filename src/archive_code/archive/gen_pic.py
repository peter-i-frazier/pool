import csv
import subprocess
def construct_stripes(img_folder, tmp_folder, lib_folder, img_list, ts, num_row, num_col, row_offset):
    for img_name in img_list:
        to_call = ['convert', '{0}/{1}.tif'.format(img_folder, img_name), '-crop', '{0}x{1}@'.format(num_col, num_row), '{0}/{1}_%d.tif'.format(tmp_folder, img_name)]
        subprocess.call(to_call)
    for row in range(num_row):
        for col in range(num_col):
            #to_call = ['convert', '{0}/ts_{1}.tif'.format(lib_folder, ts), '{0}/r_{1}.tif'.format(lib_folder, row+row_offset), '{0}/c_{1}.tif'.format(lib_folder, col)]
            to_call = ['convert']
            for img_name in img_list:
                to_call.append('{0}/{1}_{2}.tif'.format(tmp_folder, img_name, row*30+col))
            to_call.append('+append')
            to_call.append('{0}/TS_{1}_row_{2}_col_{3}.tif'.format(lib_folder, ts, row+row_offset, col))
            subprocess.call(to_call)

def reconstruct(img_folder, img_name, tmp_folder, dest_folder, num_row, num_col):
    to_call = ['convert', '{0}/{1}.tif'.format(img_folder, img_name), '-crop', '{0}x{1}@'.format(num_col, num_row), '{0}/{1}_%d.tif'.format(tmp_folder, img_name)]
    subprocess.call(to_call)
    for row in range(num_row):
        to_call = ['convert']
        for col in range(num_col):
            to_call.append('{0}/{1}_{2}.tif'.format(tmp_folder, img_name, row*30+col))
        to_call.append('+append')
        to_call.append('{0}/{1}_row_{2}.tif'.format(tmp_folder, img_name, row))
        subprocess.call(to_call)
    to_call = ['convert']
    for row in range(num_row):
        to_call.append('{0}/{1}_row_{2}.tif'.format(tmp_folder, img_name, row))
    to_call.append('-append')
    to_call.append('{0}/{1}_reconstruct.tif'.format(dest_folder, img_name))
    subprocess.call(to_call)

def rank_stripes(dest_folder, lib_folder, rank_table_folder, rank_table_name, limit=-1):
    f = open("{0}/{1}.csv".format(rank_table_folder, rank_table_name), 'rb')
    csv_f = csv.reader(f)
    for l, line in enumerate(csv_f):
        if limit > 0 and l >= limit:
            break
        if l == 1:
            to_append = "{0}/TS{1}_row_{2}_col_{3}.png".format(lib_folder, line[2], line[4], line[5])
            target = "{0}/{1}.png".format(dest_folder, rank_table_name)
            subprocess.call(['convert', to_append, '-append', target])
        elif l > 1:
            to_append = "{0}/TS{1}_row_{2}_col_{3}.png".format(lib_folder, line[2], line[4], line[5])
            subprocess.call(['convert', target, to_append, '-append', target])
    f.close()

if __name__ == "__main__":
# construct stripes
    TS1 = ['TS1_sfp_1', 'TS1_sfp_2', 'TS1_AcpS_1', 'TS1_AcpS_2']
    #TS2 = ['TS2_sfp_1', 'TS2_sfp_2']
    #TS3 = ['TS3_sfp_1', 'TS3_sfp_2', 'TS3_AcpS_1', 'TS3_AcpS_2']
    #TS4 = ['TS4_sfp_1', 'TS4_sfp_2', 'TS4_AcpS_1', 'TS4_AcpS_2']
    #TS5 = ['TS5_sfp_1', 'TS5_sfp_2', 'TS5_AcpS_1', 'TS5_AcpS_2']
    subprocess.call(['rm', '-rf', 'tmp_img'])
    subprocess.call(['mkdir', 'tmp_img'])
    construct_stripes('org_img', 'tmp_img', 'lib_img', TS1, 1, 9, 30, 1)
    #construct_stripes(org_img, "tmp_img", "lib_img", "lib_img", TS1, 1, 10, 30)
    #construct_stripes(org_img, "tmp_img", "lib_img", "lib_img", TS2, 2, 20, 30)
    #construct_stripes(org_img, "tmp_img", "lib_img", "lib_img", TS3, 3, 19, 30)
    #construct_stripes(org_img, "tmp_img", "lib_img", "lib_img", TS4, 4, 20, 30)
    #construct_stripes(org_img, "tmp_img", "lib_img", "lib_img", TS5, 5, 20, 30)
# rank stripes
    #subprocess.call(['rm', '-rf', 'result_img'])
    #subprocess.call(['mkdir', 'result_img'])
    #print "measure 1"
    #rank_stripes("result_img", "lib_img", "rank_table", "rank_by_measure_1", limit=300)
    #print "measure 2"
    #rank_stripes("result_img", "lib_img", "rank_table", "rank_by_measure_2", limit=300)
    #print "measure 3"
    #rank_stripes("result_img", "lib_img", "rank_table", "rank_by_measure_3", limit=300)
    #print "measure 4"
    #rank_stripes("result_img", "lib_img", "rank_table", "rank_by_measure_4", limit=300)

    #tablename = "TS5_sfp_1_rank"
    #print tablename
    #rank_stripes("result_img", "lib_img", "rank_table", tablename)
    #tablename = "TS5_sfp_2_rank"
    #print tablename
    #rank_stripes("result_img", "lib_img", "rank_table", tablename)
    #tablename = "TS5_AcpS_1_rank"
    #print tablename
    #rank_stripes("result_img", "lib_img", "rank_table", tablename)
    #tablename = "TS5_AcpS_2_rank"
    #print tablename
    #rank_stripes("result_img", "lib_img", "rank_table", tablename)
    #reconstruct('org_img', 'TS1_3', 'tmp_img', 'result_img', 9, 30)
