import sys
import os
def testsys(filename,filenameout,filenameout1):
    with open('inter-species_gene_synteny_file', 'r', encoding='utf-8') as fu:
        content_unique = fu.readlines()
    unique_1, unique_2, list2 = [], [], []
    for lines in content_unique:
        line = lines.rstrip('\n').split('\t')
        unique_1.append(line[0])
        unique_2.append(line[1])
        list2.extend(line)
    print(len(set(unique_1)), len(set(unique_2)), len(set(list2)), end='\t')
    with open(filename) as fm:
        content_merge = fm.readlines()
    merge = []
    for ele_con_mer in content_merge[1:]:
        ele_mer = ele_con_mer.rstrip('\n').split('\t')
        merge.append(ele_mer)
    with open('gff_File') as fg:
        content_gff = fg.readlines()
    gff = []
    for ele_con_gff in content_gff:
        ele_gff = ele_con_gff.rstrip('\n').split('\t')
        gff.append(ele_gff)
    print(len(gff))
# list：merge、gff、unique1、2

    with open(filenameout, 'a') as fp:
        for i in merge:
            big_name, big_left, big_right, typ = i[0], int(i[1]), int(i[2]), i[3]
            for j in gff:
                if big_name in j and int(j[2]) > big_left and int(j[3]) < big_right:
                    small_name, small_left, small_right = j[1], j[2], j[3]

                    if small_name.startswith('species_prefix1'):
                        if small_name in unique_1:
                            fp.write(str(big_name)+'\t'+str(big_left)+'\t'+str(big_right)+'\t')
                            fp.write(str(small_name)+'\t'+str(small_left)+'\t'+str(small_right)+'\t'+str(typ)+'\t')
                            other_name = unique_2[unique_1.index(small_name)]
                            fp.write(str(other_name)+'\t')
                            for other_gff in gff:
                                if other_name in other_gff:
                                    other_left, other_right, other_big_name = int(other_gff[2]), int(other_gff[3]), other_gff[0]
                                    fp.write(str(other_left)+'\t'+str(other_right)+'\t')
                                    for other_merge in merge:
                                        if other_big_name in other_merge and other_left < int(other_merge[2]) and other_right > int(other_merge[1]):
                                            other_type = other_merge[3]
                                            fp.write(str(other_type)+'')
                                    fp.write('\n')

                    if small_name.startswith('species_prefix2'):
                        if small_name in unique_2:
                            fp.write(str(big_name) + '\t' + str(big_left) + '\t' + str(big_right) + '\t')
                            fp.write(str(small_name) + '\t' + str(small_left) + '\t' + str(small_right)+'\t'+str(typ)+'\t')
                            other_name = unique_1[unique_2.index(small_name)]
                            fp.write(str(other_name) + '\t')
                            for other_gff in gff:
                                if other_name in other_gff:
                                    other_left, other_right, other_big_name = int(other_gff[2]), int(other_gff[3]), other_gff[0]
                                    fp.write(str(other_left) + '\t' + str(other_right) + '\t')
                                    for other_merge in merge:
                                        if other_big_name in other_merge and other_left < int(other_merge[2]) and other_right > int(other_merge[1]):
                                            other_type = other_merge[3]
                                            fp.write(str(other_type)+'')
                                    fp.write('\n')

    with open(filenameout, 'r', encoding='utf-8') as ft:
        test_content = ft.readlines()
    with open(filenameout1,'a', encoding='utf-8') as fr:
        for test_line in test_content:
            test_list = test_line.rstrip('\n').split('\t')
            if 'same' not in test_list and 'same' not in test_list[10]:
                for wr in test_list:
                    fr.write(str(wr) + '\t')
                fr.write('\n')
    print('process finished!')
if __name__ == '__main__':
    try:
        filename = sys.argv[1]
        filenameout = sys.argv[2]
        filenameout1 = sys.argv[3]
        testsys(filename,filenameout,filenameout1)
    except Exception as e:
        print(sys.argv)
        print(os.getcwd())
        print(e)
