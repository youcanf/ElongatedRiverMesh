#// writeFort14_1direction 
#// Author: Youcan Feng
#// Created on Apr 8, 2020


import os
import numpy as np
import traceback


xs_start = [1, 2, 3, 4, 5, 6]   
width_elems = 5 #width_nodes == width_elems + 1
num_XS = 229 
num_XS_estuary = 174 
flag_weir_extra_rows = True
             
def write_1dir_weir24_Fort14_customerized(in_file, out_file):
    
    width_node = width_elems + 1 
    center_node_1st = int(width_node / 2) 
    center_node_2nd = int(width_node / 2) + 1
    list_elems_new = [0, 0, 0]
    if flag_weir_extra_rows:
        temp_elems = width_elems - 2
    else:
        temp_elems = width_elems
    
    flag_odd = 0; 
    if width_elems % 2 > 0: 
        flag_odd = 1
    flag_even = 1 - flag_odd     
    
    with open(in_file, 'r') as f_in:
        f_in.readline()
        [ne, nn] = f_in.readline().strip().split()[0:2]
        ne = int(ne); nn = int(nn)
        for i in range(nn):
            f_in.readline()
        if not (ne == 1 + (num_XS-1)*temp_elems + num_XS_estuary-1):
            print('XS nums do not match fort.14.')
        else:
            for j in range((num_XS-1)*temp_elems):
                
                temp = f_in.readline()
                [id, num, node1, node2, node3, node4] = [int(temp2) for temp2 in temp.strip().split()]       
                rank_node1 = int(node1/width_node) if node1%width_node != 0 else int(node1/width_node)-1
                rank_node2 = int(node2/width_node) if node2%width_node != 0 else int(node2/width_node)-1
                rank_node3 = int(node3/width_node) if node3%width_node != 0 else int(node3/width_node)-1
                rank_node4 = int(node4/width_node) if node4%width_node != 0 else int(node4/width_node)-1

                uni_node1 = node1%width_node if node1%width_node != 0 else width_node
                uni_node2 = node2%width_node if node2%width_node != 0 else width_node
                uni_node3 = node3%width_node if node3%width_node != 0 else width_node
                uni_node4 = node4%width_node if node4%width_node != 0 else width_node
                
                test_pos = 0
                dist = 0
                if(uni_node1 < center_node_1st-dist): test_pos = test_pos + 1
                if(uni_node2 < center_node_1st-dist): test_pos = test_pos + 1
                if(uni_node3 < center_node_1st-dist): test_pos = test_pos + 1
                if(uni_node4 < center_node_1st-dist): test_pos = test_pos + 1
                
                if(uni_node1 > center_node_2nd+dist): test_pos = test_pos - 1
                if(uni_node2 > center_node_2nd+dist): test_pos = test_pos - 1
                if(uni_node3 > center_node_2nd+dist): test_pos = test_pos - 1
                if(uni_node4 > center_node_2nd+dist): test_pos = test_pos - 1            
                
                if (not flag_weir_extra_rows) and (test_pos > 0):
                    if (uni_node2 == center_node_1st and uni_node3 == center_node_1st):
                        continue  
                elif (not flag_weir_extra_rows) and (test_pos < 0):
                    if (uni_node1 == center_node_2nd and uni_node4 == center_node_2nd):
                        continue  
                else:    
                    if(rank_node1 == rank_node2):
                        if(uni_node2 == uni_node3):
                            list_elems_new = np.vstack([list_elems_new, [node1, node2, node3]])
                            list_elems_new = np.vstack([list_elems_new, [node1, node3, node4]])
                        else:
                            print('error at nodes %s %s %s %s\n' %(node1, node2, node3, node4))      
            temp = f_in.readline()
            [id, num, node1, node2, node3, node4] = [int(temp2) for temp2 in temp.strip().split()]
            list_elems_new = np.vstack([list_elems_new, [node1, node2, node3]])
            list_elems_new = np.vstack([list_elems_new, [node1, node3, node4]])
            
            for j in range(num_XS_estuary-1):
                temp = f_in.readline()
                [id, num, node1, node2, node3, node4] = [int(temp2) for temp2 in temp.strip().split()]
                list_elems_new = np.vstack([list_elems_new, [node1, node2, node3]])
                list_elems_new = np.vstack([list_elems_new, [node1, node3, node4]])
    list_elems_new = np.delete(list_elems_new, 0, 0)
    
    with open(in_file, 'r') as f_in:
        with open(out_file, 'w') as f_out:
            f_out.write(f_in.readline())
            f_in.readline()
            f_out.write('%s\t%s\n' %(len(list_elems_new), nn))

            for i in range(nn):
                f_out.write(f_in.readline())
            
            for i in range(ne):
                next(f_in) 
                
            for i in range(int(len(list_elems_new))):
                f_out.write('\t%s\t%s\t%s\t%s\t%s\n' %(i+1, 3, list_elems_new[i][0], list_elems_new[i][1], list_elems_new[i][2]))
                     
            for line in f_in.readlines():
                f_out.write(line)                            
                                             
def main():
    
    in_file = 'Outputs/fort14_quadrils_interpolated.14'
    out_file = 'Outputs/fort_3elems_infill_1dir.14'
    write_1dir_weir24_Fort14_customerized(in_file, out_file)  

if __name__ == '__main__':
    main()
