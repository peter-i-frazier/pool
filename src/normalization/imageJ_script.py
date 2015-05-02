from ij import IJ
import csv
import os

def construct_stripes(img_folder, lib_folder, img_list, ts, num_row, num_col, row_offset):
    imp_list = []
    imp_row_list = []
    imp_col_list = []
    imp_dup_list = []
    w = []
    h = []
    for k, img_name in enumerate(img_list):
        imp_list.append(IJ.openImage(img_folder + "/" + img_name + ".tif"))
        size = imp_list[k].getDimensions()
        w.append(size[0])
        h.append(size[1])
    # for r in range(num_row):
    # 	imp_row_list.append(IJ.openImage(lib_folder+"/row_"+str(r+row_offset)+".tif"))
    # for c in range(num_col):
    # 	imp_col_list.append(IJ.openImage(lib_folder+"/col_"+str(c)+".tif"))
    for r in range(num_row):
        for c in range(num_col):
            # IJ.run(imp_row_list[r], "Duplicate...", "title=row")
            # IJ.run(imp_col_list[c], "Duplicate...", "title=col")
            for k in range(len(img_list)):
            	imp_list[k].setRoi(w[k] / num_col * c, h[k] / num_row * r, w[k] / num_col, h[k] / num_row)
                IJ.run(imp_list[k], "Duplicate...", "title=" + str(k))
            IJ.run("Combine...", "stack1=0 stack2=1")
            for k in range(2, len(img_list)):
                IJ.run("Combine...", "stack1=[Combined Stacks] stack2="+str(k))
            imp = IJ.getImage()
            IJ.saveAs(imp, "TIFF", lib_folder + "/TS" + str(ts) + "_" + str(r+row_offset) + "_" + str(c) + ".tif")
            imp.close()


def parse(line):
    ts = line[0]
    row = "ABCDEFGHIJKLMNOPQRST".index(line[1][0])
    col = int(line[1][1:]) - 1
    return (ts, row, col)

def rank_stripes(table_folder, table_name, lib_folder, dest_folder):
    f = open(table_folder + "/" + table_name + ".csv", 'rb')
    csv_f = csv.reader(f)
    flag = False
    for l, line in enumerate(csv_f):
        if l == 1:
            line1 = line
        if l == 2:
            flag = True
            ts, r, c = parse(line1)
            imp = IJ.openImage(lib_folder+"/TS" + str(ts) + "_"+str(r)+"_"+str(c)+".tif")
            IJ.run(imp, "Duplicate...", "title=1")
            ts, r, c = parse(line)
            imp = IJ.openImage(lib_folder+"/TS" + str(ts) + "_"+str(r)+"_"+str(c)+".tif")
            IJ.run(imp, "Duplicate...", "title=2")
            IJ.run("Combine...", "stack1=1 stack2=2 combine")
        if l > 2:
            ts, r, c = parse(line)
            imp = IJ.openImage(lib_folder+"/TS" + str(ts) + "_"+str(r)+"_"+str(c)+".tif")
            IJ.run(imp, "Duplicate...", "title=2")
            IJ.run("Combine...", "stack1=[Combined Stacks] stack2=2 combine")
    if flag:
        imp = IJ.getImage()
        IJ.saveAs(imp, "TIFF", dest_folder + "/" + table_name + ".tif")
        imp.close()
    f.close()


root = "/Users/jialeiwang/"
org_img = root + "peptide-catalysis/src/normalization/org_img"
lib_img = root + "peptide-catalysis/src/normalization/lib_img"
test_img = root + "peptide-catalysis/src/normalization/test_img"
table_folder = root + "peptide-catalysis/src/normalization/rank_table"
dup_img = root + "peptide-catalysis/src/normalization/dup_img"
rank_img = root + "peptide-catalysis/src/normalization/rank_img"

TS1 = ['TS1_sfp_1', 'TS1_sfp_2', 'TS1_AcpS_1', 'TS1_AcpS_2']
TS2 = ['TS2_sfp_1', 'TS2_sfp_2']
TS3 = ['TS3_sfp_1', 'TS3_sfp_2', 'TS3_AcpS_1', 'TS3_AcpS_2']
TS4 = ['TS4_sfp_1', 'TS4_sfp_2', 'TS4_AcpS_1', 'TS4_AcpS_2']
TS5 = ['TS5_sfp_1', 'TS5_sfp_2', 'TS5_AcpS_1', 'TS5_AcpS_2']
#construct_stripes(org_img, lib_img, TS3, 3, 19, 30, 0)
#construct_stripes(org_img, lib_img, TS4, 4, 20, 30, 0)
#construct_stripes(org_img, lib_img, TS5, 5, 20, 30, 0)
#rank_stripes(lib_img, test_img, 3, 19, 30, 0)
#rank_stripes(lib_img, test_img, 4, 20, 30, 0)
#rank_stripes(lib_img, test_img, 5, 20, 30, 0)
#dup_idx = [1, 2, 5, 6, 107, 146, 159, 276, 277, 280, 281, 451, 453, 469, 522, 556, 577, 627, 688, 840, 841, 844, 845, 846, 847, 851, 859, 911, 972, 991, 1009, 1060, 1124, 1144, 1180, 1306, 1375, 1386, 1398, 1524, 1904]
#for i in dup_idx:
#    rank_stripes(table_folder, 'sfp_1_dup_' + str(i), lib_img, dup_img)
for filename in os.listdir(table_folder):
    if filename != '.DS_Store':
        filename = filename[0:(len(filename) - 4)]
        print filename
        rank_stripes(table_folder, filename, lib_img, rank_img)
