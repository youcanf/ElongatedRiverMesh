#// Flowlines_BuildXS_weir24_customerized
#// Author: Youcan Feng
#// Create Date: 5/7/2020
#// Purpose: prototype to build trapezoidal channel XS based on NHD streamline

import os
import numpy as np
import ogr
import gc



class NHDFlowline(object):

    #number of elements cross the channel
    num_elems = 5
    
    flag_flat = False 
    flag_flat_startend = False
    S_hypo_startend = 0.01/100  
    
    ### debug only
    ### idealized case with a constant slope
    flag_ideal = False   #master key
    BIAS_elev_m = 10 
    flag_uniform_width = False
    bank_wall = 10 #default 0 meaning turning it off
    #overwrite
    if flag_ideal:
        flag_flat_startend = True
    ###
    
    flag_weir = True #weir height must be higher than bathy
    weir_type = 24 #23 #NB: if 24, make sure channel is at least 3 elems wide
    flag_weir_extra_rows = True #there are extra elements outside the weirs
          
    flag_smoothchannel = False #True
    flag_smoothbathy = False #True
    flag_smoothwidth = False #True
    flag_smoothbank = False #True
    
    DX = 100 #80, 100 #m
    DTHETALIM = 20 #deg
    WIDTHLIM = 0  #0 #m
    H0 = 0.01 #1cm
     
    #(2)Flowlines (you probably need the file with the sorted connectivity
    flag_XS_file = 2 #2: HECRAS
    file_path = 'Outputs/NHD_XS_pts_split2_joined_sorted_half.shp'
    outfort14file = "Outputs/fort14_quadrils_interpolated.14"
   

    def __init__(self):
        self.list_ND_Geo_master=[]
        
        if self.flag_uniform_width:
            self.constantwidth = 10 #m
            
        if self.flag_smoothchannel:
            self.channelneedsmooth = True
            self.mindx = 35 #-9999      #minimum distance between two XS; the previous XS will be deleted if the distance is smaller than mindx

        if self.flag_smoothbathy:
            self.bathyneedsmooth = True
            self.mindh = 1 #m #minimum bathy difference within the moving window 
            
        if self.flag_smoothwidth:
            self.widthneedsmooth = True
            self.mindw = 0.5 #m #minimum width difference within the moving window
            
        if self.flag_smoothbank:
            self.bankneedsmooth = True
            self.mindhb = 0.5 
            
        if self.flag_weir and self.weir_type==24:
            self.weir24_extra_width = 1  
            self.weir24_elev_diff = 0.25  
            self.weir24_ending_elev = 0.25       #customerized
            self.weir24_ending_XS_flag = False    #customerized
            self.weir24_fake_botelev = -999        #customerized
            if self.flag_weir_extra_rows:
                self.weir24_extra_rows_width = 10
                        
    def CPP_xy(self, lon_deg, lat_deg):
        
        Rearth = 6378206.4
        SLAM0 = -79.0
        SFEA0 = 35.0
        
        x = Rearth*(lon_deg-SLAM0)*np.pi/180.*np.cos(SFEA0*np.pi/180.)
        y = lat_deg*np.pi/180.*Rearth
        
        return (x,y)  
    
    def CPP_lonlat(self, x, y):
        
        Rearth = 6378206.4
        SLAM0 = -79.0
        SFEA0 = 35.0  
          
        lon_deg = SLAM0 + x /np.cos(SFEA0*np.pi/180.)/Rearth/np.pi*180.
        lat_deg = y /Rearth/np.pi*180.
        
        return(lon_deg, lat_deg)
    
    def Turn_deg(self, list_xy):  
                                 
        exterior_angle_deg = np.pi
        
        x1 = list_xy[0][0]; y1 = list_xy[0][1]
        x2 = list_xy[1][0]; y2 = list_xy[1][1]
        x3 = list_xy[2][0]; y3 = list_xy[2][1]
      
        check_left = (y3-y1)*(x2-x1) - (x3-x1)*(y2-y1);
        check_left = 1 if check_left > 0 else -1
        
        check_obtuse = (((y2-y1)**2+(x2-x1)**2) + ((y3-y2)**2+(x3-x2)**2)) - ((y3-y1)**2+(x3-x1)**2)
        check_obtuse = 1 if check_obtuse < 0 else -1
          
        a2 = (y2-y1)**2+(x2-x1)**2
        b2 = (y3-y2)**2+(x3-x2)**2
        c2 = (y3-y1)**2+(x3-x1)**2
        try:
            angle_cos = (a2+b2-c2) / (2*np.sqrt(a2)*np.sqrt(b2))
            if ( (angle_cos > 1) or (angle_cos < -1) ):
                pass
            else:
                exterior_angle_deg = np.arccos( angle_cos )
        except:
            print ('Error in Turn_deg, x=%s, y=%s, a2=%s, b2=%s, c2=%s, cos=%s' %(x2,y2,a2,b2,c2,angle_cos))
        exterior_angle_deg = np.pi - exterior_angle_deg
        exterior_angle_deg = check_left * exterior_angle_deg * 180 / np.pi
        
        return exterior_angle_deg
                   
    def solver_bisect(self, data):
        
        k_perpend = data[0]
        b_perpend = data[1]
        x1 = data[2]; y1 = data[3]
        x2 = data[4]; y2 = data[5]
        x3 = data[6]; y3 = data[7]
        check_left = data[8]
        
        ymin = min(y1, y3, y2)
        ymax = max(y1, y3, y2)

        rtol = 1e-5
        maxiter = 5000
        
        it = 0; xmid = x2; ymid = y2
        errorfun = lambda x, y: self.Turn_deg([[x1,y1],[x,y],[x3,y3]]) - check_left * self.DTHETALIM
        
        residue = np.abs( errorfun(x2,y2) )
        while (residue > rtol and it < maxiter): 
            ymid = (ymin + ymax) / 2.
            if k_perpend != -9999: 
                xmin = (ymin - b_perpend) / k_perpend
                xmid = (ymid - b_perpend) / k_perpend
            else: 
                xmin = x2
                xmid = x2
            f_lower = errorfun(xmin, ymin)
            f_mid = errorfun(xmid, ymid)
            
            if (f_lower*f_mid < 0): 
                ymax = ymid
            elif (f_lower*f_mid > 0): 
                ymin = ymid
            else:
                if k_perpend != -9999: xmid = (ymid - b_perpend) / k_perpend
                else: xmid = x2
                break

            residue = np.abs( f_mid )
            if not self.channelneedsmooth:
                self.channelneedsmooth = True  
            it = it + 1
        
        return xmid, ymid 
            
    def Solve_xy(self, list_xy):
         
        k_perpend = -9999; b_perpend = 0      
        x1 = list_xy[0][0]; y1 = list_xy[0][1]
        x2 = list_xy[1][0]; y2 = list_xy[1][1]
        x3 = list_xy[2][0]; y3 = list_xy[2][1]
        
        centroid_quadril_1_x = (x1+x2)/2.; centroid_quadril_1_y = (y1+y2)/2.;
        centroid_quadril_2_x = (x2+x3)/2.; centroid_quadril_2_y = (y2+y3)/2.;
           
        Turn_deg_orig = self.Turn_deg([[x1,y1],[x2,y2],[x3,y3]])
        check_left = 0
        if Turn_deg_orig > 0: check_left = 1
        elif Turn_deg_orig < 0: check_left = -1
                       
        if (centroid_quadril_2_y - centroid_quadril_1_y) != 0:
            k_perpend = -(centroid_quadril_2_x - centroid_quadril_1_x) / (centroid_quadril_2_y - centroid_quadril_1_y)
            b_perpend = y2 - k_perpend * x2
        
        data = [k_perpend, b_perpend, x1, y1, x2, y2, x3, y3, check_left]
        new_x, new_y = self.solver_bisect(data)
        return new_x, new_y
        
    def SmoothChannel(self):
        
        list_ND_Geo_3 = self.list_ND_Geo_master
        self.channelneedsmooth = False
        rtol = 1e-1
        list_to_remove = []

        num_XS = int(len(list_ND_Geo_3))
    
        for i in range(num_XS):
            previous_lon = list_ND_Geo_3[i-1][0] if i > 0 else list_ND_Geo_3[i][0]
            previous_lat = list_ND_Geo_3[i-1][1] if i > 0 else list_ND_Geo_3[i][1]
            start_lon = list_ND_Geo_3[i][0] 
            start_lat = list_ND_Geo_3[i][1]  
            next_lon = list_ND_Geo_3[i+1][0] if i < num_XS-1 else list_ND_Geo_3[i][0]
            next_lat = list_ND_Geo_3[i+1][1] if i < num_XS-1 else list_ND_Geo_3[i][1]
            (previous_x, previous_y) = self.CPP_xy(previous_lon, previous_lat)
            (start_x, start_y) = self.CPP_xy(start_lon, start_lat)
            (next_x, next_y) = self.CPP_xy(next_lon, next_lat)  
            
            if (i > 0 and i < num_XS-1):
                list_xy = [[previous_x, previous_y], [start_x, start_y], [next_x, next_y]]
                ang_deg = self.Turn_deg(list_xy)
                if ( (ang_deg > (self.DTHETALIM + rtol)) or (ang_deg < (-self.DTHETALIM - rtol)) ):
                    start_x, start_y = self.Solve_xy(list_xy)
                    start_lon, start_lat = self.CPP_lonlat(start_x, start_y)
                    list_ND_Geo_3[i][0] = start_lon 
                    list_ND_Geo_3[i][1] = start_lat 
                    if ( np.sqrt((start_x - next_x)**2 + (start_y - next_y)**2) < self.mindx):
                        list_to_remove.append(i)
                     
        for i in sorted(list_to_remove, reverse=True):
            del list_ND_Geo_3[i]
            print ('XS=%s is removed.' %(i))
        
        self.list_ND_Geo_master = list_ND_Geo_3            
     
     
    def SmoothBathy(self):
        
        self.bathyneedsmooth = False
        list_ND_Geo_3 = self.list_ND_Geo_master
        num_XS = int(len(list_ND_Geo_3))
        kernal = [-42, -21, -2, 15, 30, 43, 54, 63, 70, 75, 78, 79, 78, 75, 70, 63, 54, 43, 30, 15, -2, -21, -42]        
        half_win = int(len(kernal)/2.)
        sum_k = sum(kernal)
        max_diff = 0        
        for i in range(half_win, num_XS-int(half_win+1)):
            list_bathy = [list_ND_Geo_3[j][2] for j in range(i-half_win, i+int(half_win+1))]
            if min(list_bathy) != self.weir24_fake_botelev: #customerized
                max_diff = max(max_diff, abs(max(list_bathy) - min(list_bathy)))
                list_ND_Geo_3[i][2] = np.inner(list_bathy, kernal) / sum_k
        self.list_ND_Geo_master = list_ND_Geo_3
        if not self.bathyneedsmooth:
            if max_diff > self.mindh:
                self.bathyneedsmooth = True
                
        return max_diff
                
                
    def SmoothWidth(self):
        
        self.widthneedsmooth = False
        list_ND_Geo_3 = self.list_ND_Geo_master
        num_XS = int(len(list_ND_Geo_3))

        kernal = [-42, -21, -2, 15, 30, 43, 54, 63, 70, 75, 78, 79, 78, 75, 70, 63, 54, 43, 30, 15, -2, -21, -42]
        half_win = int(len(kernal)/2.)
        sum_k = sum(kernal)
        max_diff = 0       
        for i in range(half_win, num_XS-int(half_win+1)):
            list_width = [list_ND_Geo_3[j][5] for j in range(i-half_win, i+int(half_win+1))]
            max_diff = max(max_diff, abs(max(list_width) - min(list_width)))
            list_ND_Geo_3[i][5] = np.inner(list_width, kernal) / sum_k
            list_width = [list_ND_Geo_3[j][6] for j in range(i-half_win, i+int(half_win+1))]
            max_diff = max(max_diff, abs(max(list_width) - min(list_width)))
            list_ND_Geo_3[i][6] = np.inner(list_width, kernal) / sum_k
            
        self.list_ND_Geo_master = list_ND_Geo_3
        if not self.widthneedsmooth:
            if max_diff > self.mindw:
                self.widthneedsmooth = True
                                
        return max_diff
    

    def SmoothBank(self):
        
        self.bankneedsmooth = False
        list_ND_Geo_3 = self.list_ND_Geo_master
        num_XS = int(len(list_ND_Geo_3))
        
        kernal = [-42, -21, -2, 15, 30, 43, 54, 63, 70, 75, 78, 79, 78, 75, 70, 63, 54, 43, 30, 15, -2, -21, -42]
        half_win = int(len(kernal)/2.)
        sum_k = sum(kernal)
        max_diff = 0     
        for i in range(half_win, num_XS-int(half_win+1)):
            list_bank = [list_ND_Geo_3[j][3] for j in range(i-half_win, i+int(half_win+1))]
            if min(list_bank) != self.weir24_fake_botelev: #customerized
                max_diff = max(max_diff, abs(max(list_bank) - min(list_bank)))
                list_ND_Geo_3[i][3] = np.inner(list_bank, kernal) / sum_k
            list_bank = [list_ND_Geo_3[j][4] for j in range(i-half_win, i+int(half_win+1))]
            if min(list_bank) != self.weir24_fake_botelev: #customerized
                max_diff = max(max_diff, abs(max(list_bank) - min(list_bank)))
                list_ND_Geo_3[i][4] = np.inner(list_bank, kernal) / sum_k
            
        self.list_ND_Geo_master = list_ND_Geo_3
        if not self.bankneedsmooth:
            if max_diff > self.mindhb:
                self.bankneedsmooth = True
                                
        return max_diff
    
                       
    def main(self):    
        
        S_hypo = 0
        bot_elev_m_first = -9999;
        bank_left_elev_m_first = -9999;
        bank_right_elev_m_first = -9999;
        left_width_m_avg = -9999;
        right_width_m_avg = -9999;
     
        driver = ogr.GetDriverByName('ESRI Shapefile')
        if self.flag_XS_file == 2:
            dataSource = driver.Open(self.file_path,0)
        layer = dataSource.GetLayer()
        num_reaches = layer.GetFeatureCount()
        spatialRef = layer.GetSpatialRef()
        list_landleftBC=[]; list_landrightBC=[]; list_channelND=[]
        list_VTC_Lon=[]; list_VTC_Lat=[]; list_anchor=[]
        list_ND_Geo=[]; list_ND_Length=[]
        list_alt_raw=[]; list_M_calib=[]
        
        #(-1) Get the total length for idealized case
        #########################################################################
        if self.flag_ideal:
            count = 0; dist_sum = 0
            bot_elev_m_last = -9999;
            bank_left_elev_m_last = -9999;
            bank_right_elev_m_last = -9999;
            left_width_m_last = -9999;
            right_width_m_last = -9999;
            left_width_m_first = -9999;
            right_width_m_first = -9999;
            
            for inFeature in layer:
                if (count + 1 == num_reaches):
                    bot_elev_m_last = inFeature.GetField('Bot_elev_m')
                    bank_left_elev_m_last = inFeature.GetField('Bank_left_')
                    bank_right_elev_m_last = inFeature.GetField('Bank_right') 
                    left_width_m_last = inFeature.GetField('Left_width')
                    right_width_m_last = inFeature.GetField('Right_widt')          
                elif (count == 0):
                    left_width_m_first = inFeature.GetField('Left_width')
                    right_width_m_first = inFeature.GetField('Right_widt') 
                                    
                lines = inFeature.GetGeometryRef()
                flag_new_feature = 1
        
                for j in range(lines.GetPointCount()-1):
                    (start_x, start_y) = self.CPP_xy(lines.GetX(j), lines.GetY(j))
                    (next_x, next_y) = self.CPP_xy(lines.GetX(j+1), lines.GetY(j+1))
                    dist = np.sqrt((next_x-start_x)**2. + (next_y-start_y)**2.)         
                    dist_sum = dist_sum + dist
                count = count + 1   
            
            S_hypo = -self.BIAS_elev_m / dist_sum
            bot_elev_m_first = bot_elev_m_last + self.BIAS_elev_m;
            if self.bank_wall > 0:
                bank_left_elev_m_first = bot_elev_m_first + self.bank_wall
                bank_right_elev_m_first = bot_elev_m_first + self.bank_wall
            else:
                bank_left_elev_m_first = bank_left_elev_m_last + self.BIAS_elev_m;
                bank_right_elev_m_first = bank_right_elev_m_last + self.BIAS_elev_m;
                
            if self.flag_uniform_width:
                left_width_m_avg = self.constantwidth / 2.
                right_width_m_avg = self.constantwidth / 2.
                self.mindx = left_width_m_avg + right_width_m_avg   
 
            
        #(1) Get the list of geometry real quick
        #########################################################################
        layer.ResetReading()
        count = 0; dist_sum = 0; dist_seg_total = 0
        if self.flag_ideal:
            
            for inFeature in layer:
            
                left_width_m = inFeature.GetField('Left_width')
                right_width_m = inFeature.GetField('Right_widt')
                if count == 0:
                    bot_elev_m_old = bot_elev_m_first
                    if self.bank_wall > 0:
                        bank_left_elev_m_old = bot_elev_m_old + self.bank_wall
                        bank_right_elev_m_old = bot_elev_m_old + self.bank_wall
                    else:
                        bank_left_elev_m_old = bank_left_elev_m_first
                        bank_right_elev_m_old = bank_right_elev_m_first
                    if self.flag_uniform_width:
                        left_width_m_old = left_width_m_avg
                        right_width_m_old = right_width_m_avg
                    else:
                        left_width_m_old = left_width_m
                        right_width_m_old = right_width_m            
                
                lines = inFeature.GetGeometryRef()
                flag_new_feature = 1
        
                for j in range(lines.GetPointCount()-1):
                    (start_x, start_y) = self.CPP_xy(lines.GetX(j), lines.GetY(j))
                    (next_x, next_y) = self.CPP_xy(lines.GetX(j+1), lines.GetY(j+1))     
                    dist = np.sqrt((next_x-start_x)**2. + (next_y-start_y)**2.)                 
                    
                    list_VTC_Lon.append(lines.GetX(j))
                    list_VTC_Lat.append(lines.GetY(j))            
                    if (flag_new_feature == 1) and (j == 0): 
                        list_anchor.append(1)

                        list_ND_Geo.append([bot_elev_m_old, bank_left_elev_m_old, bank_right_elev_m_old, left_width_m_old, right_width_m_old])
                        if count > 0: 
                            list_ND_Length.append(dist_sum)
                            if (count == 1) and (not self.flag_flat_startend):
                                list_ND_Geo[0][0] = list_ND_Geo[0][0]+self.S_hypo_startend*dist_sum
                                list_ND_Geo[0][1] = list_ND_Geo[0][1]+self.S_hypo_startend*dist_sum
                                list_ND_Geo[0][2] = list_ND_Geo[0][2]+self.S_hypo_startend*dist_sum
        
                        dist_sum = 0
                        flag_new_feature = 0
                    else: list_anchor.append(0) 
                    dist_sum = dist_sum + dist
                
                dist_seg_total = dist_seg_total + dist_sum    
                bot_elev_m_old = bot_elev_m_first + S_hypo * dist_seg_total
                if self.bank_wall > 0:
                    bank_left_elev_m_old = bot_elev_m_old + self.bank_wall
                    bank_right_elev_m_old = bot_elev_m_old + self.bank_wall
                else:
                    bank_left_elev_m_old = bank_left_elev_m_first + S_hypo * dist_seg_total
                    bank_right_elev_m_old = bank_right_elev_m_first + S_hypo * dist_seg_total
                if self.flag_uniform_width:
                    left_width_m_old = left_width_m_avg
                    right_width_m_old = right_width_m_avg
                else:
                    left_width_m_old = left_width_m
                    right_width_m_old = right_width_m
                
                if (count + 1 == num_reaches): 
                    if self.flag_flat_startend:
                        list_ND_Geo.append([bot_elev_m_old, bank_left_elev_m_old, bank_right_elev_m_old, left_width_m_old, right_width_m_old])
                    else:
                        list_ND_Geo.append([bot_elev_m_old-self.S_hypo_startend*dist_sum, bank_left_elev_m_old-self.S_hypo_startend*dist_sum, bank_right_elev_m_old-self.S_hypo_startend*dist_sum, left_width_m_old, right_width_m_old])
                    list_ND_Length.append(dist_sum)  
                
                count = count + 1                
                
        else: 
            
            for inFeature in layer:
                
                if self.flag_XS_file == 2: 
                    bot_elev_m = inFeature.GetField('Bot_elev_m')
                    bank_left_elev_m = inFeature.GetField('Bank_left_')
                    bank_right_elev_m = inFeature.GetField('Bank_right')
                    left_width_m = inFeature.GetField('Left_width')
                    right_width_m = inFeature.GetField('Right_widt')
 
                if count == 0:
                    bot_elev_m_old = bot_elev_m
                    bank_left_elev_m_old = bank_left_elev_m
                    bank_right_elev_m_old = bank_right_elev_m
                    left_width_m_old = left_width_m
                    right_width_m_old = right_width_m            
                
                lines = inFeature.GetGeometryRef()
                flag_new_feature = 1
        
                for j in range(lines.GetPointCount()-1):
                    (start_x, start_y) = self.CPP_xy(lines.GetX(j), lines.GetY(j))
                    (next_x, next_y) = self.CPP_xy(lines.GetX(j+1), lines.GetY(j+1))    
                    dist = np.sqrt((next_x-start_x)**2. + (next_y-start_y)**2.)
                     
                    list_VTC_Lon.append(lines.GetX(j))
                    list_VTC_Lat.append(lines.GetY(j))            
                    if (flag_new_feature == 1) and (j == 0): 
                        list_anchor.append(1)
                        list_ND_Geo.append([bot_elev_m_old, bank_left_elev_m_old, bank_right_elev_m_old, left_width_m_old, right_width_m_old])
                        if count > 0: 
                            list_ND_Length.append(dist_sum)
                            if (count == 1) and (not self.flag_flat_startend):
                                list_ND_Geo[0][0] = list_ND_Geo[0][0]+self.S_hypo_startend*dist_sum
                                list_ND_Geo[0][1] = list_ND_Geo[0][1]+self.S_hypo_startend*dist_sum
                                list_ND_Geo[0][2] = list_ND_Geo[0][2]+self.S_hypo_startend*dist_sum
                        dist_sum = 0
                        flag_new_feature = 0
                    else: list_anchor.append(0) 
                    dist_sum = dist_sum + dist
                    
                bot_elev_m_old = bot_elev_m
                bank_left_elev_m_old = bank_left_elev_m
                bank_right_elev_m_old = bank_right_elev_m
                left_width_m_old = left_width_m
                right_width_m_old = right_width_m
                
                
                if (count + 1 == num_reaches): 
                    if self.flag_flat_startend:
                        list_ND_Geo.append([bot_elev_m_old, bank_left_elev_m_old, bank_right_elev_m_old, left_width_m_old, right_width_m_old])
                    else:
                        list_ND_Geo.append([bot_elev_m_old-self.S_hypo_startend*dist_sum, bank_left_elev_m_old-self.S_hypo_startend*dist_sum, bank_right_elev_m_old-self.S_hypo_startend*dist_sum, left_width_m_old, right_width_m_old])
                    list_ND_Length.append(dist_sum)             
                
                count = count + 1    
    
                    
        #(2) Resample the centerline 
        #########################################################################
        outDriver = ogr.GetDriverByName("ESRI Shapefile")
        if os.path.exists(self.outShapefile_center):
            outDriver.DeleteDataSource(self.outShapefile_center)
        outDataSource = outDriver.CreateDataSource(self.outShapefile_center)
        outLayer = outDataSource.CreateLayer("XS", spatialRef, geom_type = ogr.wkbPoint)   
        outLayerDefn = outLayer.GetLayerDefn()    
        newField = ogr.FieldDefn('XS_ID', ogr.OFTInteger64)
        outLayer.CreateField(newField)  
        start_lon=list_VTC_Lon[0]; start_lat=list_VTC_Lat[0]
        dist_sum = 0; count = 0

        point1 = ogr.Geometry(ogr.wkbPoint)
        point1.AddPoint(start_lon, start_lat)
        outFeature = ogr.Feature(outLayerDefn)
        outFeature.SetGeometry(point1)
        outFeature.SetField('XS_ID', count)
        outLayer.CreateFeature(outFeature)
        outFeature = None     
        
        count_ND = 0; dist_seg_total = 0
        list_ND_Geo_2=[]; 
        list_ND_Geo_2.append([start_lon, start_lat, list_ND_Geo[0][0], list_ND_Geo[0][1], list_ND_Geo[0][2], list_ND_Geo[0][3], list_ND_Geo[0][4]])
        
        for i in range(len(list_VTC_Lon)):
            if list_anchor[i] == 1:
                bot_elev_m_start = list_ND_Geo[count_ND][0]
                bank_left_elev_m_start = list_ND_Geo[count_ND][1]
                bank_right_elev_m_start = list_ND_Geo[count_ND][2]
                left_width_m_start = list_ND_Geo[count_ND][3]
                right_width_m_start = list_ND_Geo[count_ND][4]
                length = list_ND_Length[count_ND]
                
                bot_elev_m_next = list_ND_Geo[count_ND+1][0]
                bank_left_elev_m_next = list_ND_Geo[count_ND+1][1]
                bank_right_elev_m_next = list_ND_Geo[count_ND+1][2]
                left_width_m_next = list_ND_Geo[count_ND+1][3]
                right_width_m_next = list_ND_Geo[count_ND+1][4]
                
                So_bot_elev = (bot_elev_m_next - bot_elev_m_start) / length
                So_bank_left_elev = (bank_left_elev_m_next - bank_left_elev_m_start) / length
                So_bank_right_elev = (bank_right_elev_m_next - bank_right_elev_m_start) / length
                So_left_width = (left_width_m_next - left_width_m_start) / length
                So_right_width = (right_width_m_next - right_width_m_start) / length
                
                count_ND = count_ND + 1
                dist_seg_total = 0 
                
            (current_vtc_lon, current_vtc_lat) = (list_VTC_Lon[i], list_VTC_Lat[i])
            (start_x, start_y) = self.CPP_xy(start_lon, start_lat)
            (current_vtc_x, current_vtc_y) = self.CPP_xy(current_vtc_lon, current_vtc_lat)          
            dist = np.sqrt((current_vtc_x-start_x)**2. + (current_vtc_y-start_y)**2.) 
           
            if (dist_sum + dist < self.DX):
                dist_sum = dist_sum + dist
                if list_anchor[i] == 1:  
                    dist_seg_total = 0
                else:
                    dist_seg_total = dist_seg_total + dist
                (start_lon, start_lat) = (current_vtc_lon, current_vtc_lat)
                (start_x, start_y) = self.CPP_xy(start_lon, start_lat)
                continue 
            while (dist_sum + dist >= self.DX): 
                lag_m = self.DX - dist_sum
                tangent = np.arctan2((current_vtc_y-start_y), (current_vtc_x-start_x)) 
                start_x = start_x + lag_m * np.cos(tangent)
                start_y = start_y + lag_m * np.sin(tangent)
                (start_lon, start_lat) = self.CPP_lonlat(start_x, start_y)
    
                point1 = ogr.Geometry(ogr.wkbPoint)
                point1.AddPoint(start_lon, start_lat)
                outFeature = ogr.Feature(outLayerDefn)
                outFeature.SetGeometry(point1)
                outFeature.SetField('XS_ID', count)
                outLayer.CreateFeature(outFeature)
                outFeature = None
                
                if list_anchor[i] == 1:  
                    dist_seg_total = 0
                else:
                    dist_seg_total = dist_seg_total + lag_m
                if dist_seg_total > length:
                    print('Error: dist_seg_total > length of the segments from the anchor point #%s.' %(i))
                
                #customerized
                if (count_ND == num_reaches): 
                    center_elev = self.weir24_fake_botelev
                    bank_left_elev = self.weir24_fake_botelev
                    bank_right_elev = self.weir24_fake_botelev
                else:
                    center_elev = bot_elev_m_start + dist_seg_total * So_bot_elev
                    bank_left_elev = bank_left_elev_m_start + dist_seg_total * So_bank_left_elev
                    bank_right_elev = bank_right_elev_m_start + dist_seg_total * So_bank_right_elev
                
                left_width = left_width_m_start + dist_seg_total * So_left_width
                right_width = right_width_m_start + dist_seg_total * So_right_width
                list_ND_Geo_2.append([start_lon, start_lat, center_elev, bank_left_elev, bank_right_elev, left_width, right_width])
                
                dist = np.sqrt((current_vtc_x-start_x)**2. + (current_vtc_y-start_y)**2.) 
                count = count + 1

                if (dist < self.DX):
                    dist_sum = dist
                    if list_anchor[i] == 1:  
                        dist_seg_total = 0
                    else:
                        dist_seg_total = dist_seg_total + dist
                    (start_lon, start_lat) = (current_vtc_lon, current_vtc_lat)
                    break
                else:
                    dist_sum = 0 
                
        outLayer=None; outLayerDefn=None; outDataSource=None; outDriver=None      
        del list_ND_Geo; gc.collect() 
        self.list_ND_Geo_master = list_ND_Geo_2
        
        #(3) Smoothing centerlines 
        #########################################################################   
        if self.flag_smoothchannel:
            it = 0
            maxiter = 200
            
            while (self.channelneedsmooth and (it < maxiter)):   
                self.SmoothChannel()
                it = it + 1
                print ('%s round of channel smoothing is completed.' %(it))
            
            list_ND_Geo_2 = self.list_ND_Geo_master
        
        
        #(4) Smoothing bathy
        #########################################################################
        if self.flag_smoothbathy:
            it = 0
            maxiter = 200
            
            while (self.bathyneedsmooth and (it < maxiter)):
                max_diff = self.SmoothBathy()
                it = it + 1
                print ('%s round of bathy smoothing is completed (max_diff=%s).' %(it, max_diff))
                
            list_ND_Geo_2 = self.list_ND_Geo_master
            
            
        #(5) Smoothing width
        #########################################################################
        if self.flag_smoothwidth:
            it = 0
            maxiter = 200
            
            while (self.widthneedsmooth and (it < maxiter)):
                max_diff = self.SmoothWidth()
                it = it + 1
                print ('%s round of width smoothing is completed (max_diff=%s).' %(it, max_diff))
                
            list_ND_Geo_2 = self.list_ND_Geo_master
        
        
        #(6) Smoothing bank
        #########################################################################
        if self.flag_smoothbank:
            it = 0
            maxiter = 200
            
            while (self.bankneedsmooth and (it < maxiter)):
                max_diff = self.SmoothBank()
                it = it + 1
                print ('%s round of bank smoothing is completed (max_diff=%s).' %(it, max_diff))
                
            list_ND_Geo_2 = self.list_ND_Geo_master        
        
        ##################################################################################################################
        ##################################################################################################################      
        
        #(7) Calculate bank locations based on channel widths 
        #########################################################################
        spatialRef = layer.GetSpatialRef()
        outDriver = ogr.GetDriverByName("ESRI Shapefile")
        if os.path.exists(self.outShapefile_bank):
            outDriver.DeleteDataSource(self.outShapefile_bank)
        outDataSource = outDriver.CreateDataSource(self.outShapefile_bank)
        outLayer = outDataSource.CreateLayer("XS", spatialRef, geom_type = ogr.wkbPoint)   
        newField = ogr.FieldDefn('XS_ID', ogr.OFTString)  
        outLayer.CreateField(newField)
        newField = ogr.FieldDefn('Elev_m', ogr.OFTReal)  
        outLayer.CreateField(newField)
        outLayerDefn = outLayer.GetLayerDefn()
        node_count = 0
        list_nodes = [0, 0, 0]
        list_XS = [0] * (self.num_elems+1) 
        list_XS_estuary = [0] * 2   
        num_XS = int(len(list_ND_Geo_2))
        list_left_weirBC = [0] * 2
        list_right_weirBC = [0] * 2
    
        for i in range(num_XS):
            previous_lon = list_ND_Geo_2[i-1][0] if i > 0 else list_ND_Geo_2[i][0]
            previous_lat = list_ND_Geo_2[i-1][1] if i > 0 else list_ND_Geo_2[i][1]
            start_lon = list_ND_Geo_2[i][0] 
            start_lat = list_ND_Geo_2[i][1] 
            next_lon = list_ND_Geo_2[i+1][0] if i < num_XS-1 else list_ND_Geo_2[i][0]
            next_lat = list_ND_Geo_2[i+1][1] if i < num_XS-1 else list_ND_Geo_2[i][1]
            (previous_x, previous_y) = self.CPP_xy(previous_lon, previous_lat)
            (start_x, start_y) = self.CPP_xy(start_lon, start_lat)
            (next_x, next_y) = self.CPP_xy(next_lon, next_lat)  
            
            XS_ID =  i
            center_elev = list_ND_Geo_2[i][2]
            bank_left_elev = list_ND_Geo_2[i][3]
            bank_right_elev = list_ND_Geo_2[i][4]
            left_width = list_ND_Geo_2[i][5]
            right_width = list_ND_Geo_2[i][6]
            perpend = np.arctan2((next_x - previous_x), -(next_y - previous_y))
            ##################################################################################################################
            list_channelND_local = []
            
            if (self.flag_weir and self.weir_type == 24):
                if self.flag_weir_extra_rows:
                    extra_left_x = start_x + (left_width+self.weir24_extra_width+left_width+right_width) * np.cos(perpend)
                    extra_left_y = start_y + (left_width+self.weir24_extra_width+left_width+right_width) * np.sin(perpend)                            
                    extra_right_x = start_x - (right_width+self.weir24_extra_width+left_width+right_width) * np.cos(perpend)
                    extra_right_y = start_y - (right_width+self.weir24_extra_width+left_width+right_width) * np.sin(perpend)
                    (extra_left_lon, extra_left_lat) = self.CPP_lonlat(extra_left_x, extra_left_y)
                    (extra_right_lon, extra_right_lat) = self.CPP_lonlat(extra_right_x, extra_right_y)
                left_top_x = start_x + (left_width+self.weir24_extra_width) * np.cos(perpend)
                left_top_y = start_y + (left_width+self.weir24_extra_width) * np.sin(perpend)                            
                right_top_x = start_x - (right_width+self.weir24_extra_width) * np.cos(perpend)
                right_top_y = start_y - (right_width+self.weir24_extra_width) * np.sin(perpend)
            else:
                left_top_x = start_x + left_width * np.cos(perpend)
                left_top_y = start_y + left_width * np.sin(perpend)                     
                right_top_x = start_x - right_width * np.cos(perpend)
                right_top_y = start_y - right_width * np.sin(perpend)
            (left_top_lon, left_top_lat) = self.CPP_lonlat(left_top_x, left_top_y)    
            (right_top_lon, right_top_lat) = self.CPP_lonlat(right_top_x, right_top_y)
            
            if self.flag_flat: 
                left_top_elev = center_elev
                right_top_elev = center_elev
            else:
                if self.num_elems>1 and self.flag_weir and self.weir_type == 24: 
                    left_top_elev = bank_left_elev  
                    right_top_elev = bank_right_elev 
                else:
                    left_top_elev = bank_left_elev
                    right_top_elev = bank_right_elev
            center_new_x = start_x
            center_new_y = start_y        
            
            if self.flag_XS_file == 2: 
                left_btm_x = center_new_x + left_width * np.cos(perpend)
                left_btm_y = center_new_y + left_width * np.sin(perpend)
                right_btm_x = center_new_x - right_width * np.cos(perpend)
                right_btm_y = center_new_y - right_width * np.sin(perpend)
                
            (left_btm_lon, left_btm_lat) = self.CPP_lonlat(left_btm_x, left_btm_y)
            left_btm_elev = center_elev
            (right_btm_lon, right_btm_lat) = self.CPP_lonlat(right_btm_x, right_btm_y)
            right_btm_elev = center_elev
            
            point1 = ogr.Geometry(ogr.wkbPoint)
            point1.AddPoint(left_top_lon, left_top_lat)
            outFeature = ogr.Feature(outLayerDefn)
            outFeature.SetGeometry(point1)
            outFeature.SetField('XS_ID', XS_ID)
            outFeature.SetField('Elev_m', left_top_elev)
            outLayer.CreateFeature(outFeature)
            outFeature = None
            
            point2 = ogr.Geometry(ogr.wkbPoint)
            point2.AddPoint(left_btm_lon, left_btm_lat)
            outFeature = ogr.Feature(outLayerDefn)            
            outFeature.SetGeometry(point2)
            outFeature.SetField('XS_ID', XS_ID)
            outFeature.SetField('Elev_m', left_btm_elev) 
            outLayer.CreateFeature(outFeature)
            outFeature = None           
            
            point3 = ogr.Geometry(ogr.wkbPoint)
            point3.AddPoint(right_btm_lon, right_btm_lat)
            outFeature = ogr.Feature(outLayerDefn)  
            outFeature.SetGeometry(point3)
            outFeature.SetField('XS_ID', XS_ID)
            outFeature.SetField('Elev_m', right_btm_elev)  
            # Add created feature into layer
            outLayer.CreateFeature(outFeature)            
            outFeature = None           
            
            point4 = ogr.Geometry(ogr.wkbPoint)
            point4.AddPoint(right_top_lon, right_top_lat)
            outFeature = ogr.Feature(outLayerDefn)              
            outFeature.SetGeometry(point4)
            outFeature.SetField('XS_ID', XS_ID)
            outFeature.SetField('Elev_m', right_top_elev)  
            outLayer.CreateFeature(outFeature)            
            outFeature = None   
            
            #customerized
            if (bank_left_elev <= self.weir24_ending_elev or bank_right_elev <= self.weir24_ending_elev):
                self.weir24_ending_XS_flag = True 
                    
            # summarize for fort.14
            if (self.num_elems > 1 and self.flag_weir and self.weir_type==24):
                if not self.weir24_ending_XS_flag:
                    if self.flag_weir_extra_rows:
                        list_nodes = np.vstack([list_nodes, [extra_left_lon, extra_left_lat, left_top_elev]])
                        node_count = node_count + 1; node_extra_left = node_count
                    list_nodes = np.vstack([list_nodes, [left_top_lon, left_top_lat, left_top_elev]])
                    node_count = node_count + 1; node1 = node_count
                    list_channelND_local.append(node1) 
            list_nodes = np.vstack([list_nodes, [left_btm_lon, left_btm_lat, left_btm_elev]])
            node_count = node_count + 1; node2 = node_count
            list_channelND_local.append(node2)
            list_nodes = np.vstack([list_nodes, [right_btm_lon, right_btm_lat, right_btm_elev]])
            node_count = node_count + 1; node3 = node_count
            list_channelND_local.append(node3)
            if (self.num_elems > 1 and self.flag_weir and self.weir_type==24):
                if not self.weir24_ending_XS_flag:                   
                    list_nodes = np.vstack([list_nodes, [right_top_lon, right_top_lat, right_top_elev]])
                    node_count = node_count + 1; node4 = node_count
                    list_channelND_local.append(node4)
                    if self.flag_weir_extra_rows:
                        list_nodes = np.vstack([list_nodes, [extra_right_lon, extra_right_lat, right_top_elev]])
                        node_count = node_count + 1; node_extra_right = node_count        
            if (self.num_elems > 1 and self.flag_weir and self.weir_type==24):
                if not self.weir24_ending_XS_flag: 
                    if self.flag_weir_extra_rows:
                        list_XS = np.vstack([list_XS, [node_extra_left, node1, node2, node3, node4, node_extra_right]])
                    else:
                        list_XS = np.vstack([list_XS, [node1, node2, node3, node4]])
                else:
                    list_XS_estuary = np.vstack([list_XS_estuary, [node2, node3]])
            else:
                list_XS = np.vstack([list_XS, [node2, node3]])
            list_channelND.extend(list_channelND_local) 

            if (self.num_elems > 1):
                if (self.flag_weir and self.weir_type==24):
                    pass
                else:
                    list_landleftBC.append(node1)
                    list_landrightBC.append(node4)
            elif (self.num_elems == 1):
                if not self.flag_weir: 
                    list_landleftBC.append(node2)
                    list_landrightBC.append(node3)
                
            if self.flag_weir: 
                if not self.weir24_ending_XS_flag: #customerized
                    list_left_weirBC = np.vstack([list_left_weirBC, [node2, bank_left_elev+10*self.H0]])   #+10*H0    #V9
                if not self.weir24_ending_XS_flag: #customerized
                    list_right_weirBC = np.vstack([list_right_weirBC, [node3, bank_right_elev+10*self.H0]])    #+10*H0    #V9
            print ('finish ... %s ... XS' %(i))
        outLayer = None;

        list_nodes = np.delete(list_nodes, 0, 0)
        list_XS = np.delete(list_XS, 0, 0)
        list_XS_estuary = np.delete(list_XS_estuary, 0, 0)
        list_left_weirBC = np.delete(list_left_weirBC, 0, 0)
        list_right_weirBC =  np.delete(list_right_weirBC, 0, 0)        
        
        
        #write fort.14
        #########################################################################
        with open(self.outfort14file, 'w') as f_out:
            if self.flag_weir_extra_rows:
                ne = int((len(list_XS)-1) * (self.num_elems-2) + len(list_XS_estuary)-1 +1)
            else:
                ne = int((len(list_XS)-1) * (self.num_elems) + len(list_XS_estuary)-1 +1)  
            f_out.write('New main\n')
            f_out.write('%s\t%s\n' %(ne, int(node_count)))
    
            for i in range(int(node_count)):
                f_out.write('%s\t%s\t%s\t%s\n' %(i+1, list_nodes[i][0], list_nodes[i][1], -list_nodes[i][2]))
            
            count_elems = 1
            for i in range(len(list_XS)-1):
                for j in range(len(list_XS[i])-1):
                    if self.flag_weir_extra_rows and (j in [1, 3]):
                        continue
                    else:
                        f_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(count_elems, 4, list_XS[i][j], list_XS[i][j+1], list_XS[i+1][j+1], list_XS[i+1][j]))
                        count_elems = count_elems + 1

            if self.flag_weir_extra_rows:
                f_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(count_elems, 4, list_XS[-1][2], list_XS[-1][3], list_XS_estuary[0][1], list_XS_estuary[0][0]))
            else:
                f_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(count_elems, 4, list_XS[-1][1], list_XS[-1][2], list_XS_estuary[0][1], list_XS_estuary[0][0]))
            count_elems = count_elems + 1       
            for i in range(len(list_XS_estuary)-1):
                for j in range(len(list_XS_estuary[i])-1):
                    f_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(count_elems, 4, list_XS_estuary[i][j], list_XS_estuary[i][j+1], list_XS_estuary[i+1][j+1], list_XS_estuary[i+1][j]))
                    count_elems = count_elems + 1
            print('List_XS ends at %s, List_XS_estuary ends at %s' %(len(list_XS),len(list_XS_estuary)))
                      
            f_out.write('%s = Number of open boundaries\n' %(1)) 
            f_out.write('%s = Total number of open boundary nodes\n' %(2))
            f_out.write('%s = Number of nodes for open boundary 1\n' %(2))
            if self.num_elems>1 and self.flag_weir and self.weir_type == 24:
                if self.weir24_ending_XS_flag: #customerized
                    f_out.write('%s\n%s\n' %(list_channelND[-1], list_channelND[-2]))
                else:
                    f_out.write('%s\n%s\n' %(list_channelND[-2], list_channelND[-3]))
            else:
                f_out.write('%s\n%s\n' %(list_channelND[-1], list_channelND[-2]))

            if not self.flag_weir:
                f_out.write('%s = Number of land boundaries\n' %(3))
                f_out.write('%s = Total number of land boundary nodes\n' %(int(len(list_landleftBC)+len(list_landrightBC)+2)))
                f_out.write('%s\t%s\n' %(2, 22))
                f_out.write('%s\n%s\n' %(list_channelND[0], list_channelND[1]))
            else: 
                f_out.write('%s = Number of land boundaries\n' %(3))
                f_out.write('%s = Total number of land boundary nodes\n' %(len(list_left_weirBC)+len(list_right_weirBC)+2))                
    
                f_out.write('%s\t%s\n' %(2, 22))
                if self.num_elems>1 and self.flag_weir and self.weir_type == 24:
                    f_out.write('%s\n%s\n' %(list_channelND[1], list_channelND[2]))
                else:
                    f_out.write('%s\n%s\n' %(list_channelND[0], list_channelND[1]))

            if not self.flag_weir:
                if (self.num_elems > 1):
                    f_out.write('%s\t%s\n' %(int(len(list_landrightBC)+2), 20))
                    f_out.write('%s\n' %(list_channelND[0]))
                else:
                    f_out.write('%s\t%s\n' %(int(len(list_landrightBC)), 20))
                for i in list_landrightBC:
                    f_out.write('%s\n' %(i))
            else:         
                if self.weir_type == 23:
                    BARLANCFSP = 1.0    
                    f_out.write('%s\t%s\n' %(int(len(list_right_weirBC)), 23))
                    for i in list_right_weirBC:
                        f_out.write('%s\t%s\t%s\n' %(int(i[0]), i[1], BARLANCFSP))  #list_nodes[i-1][2] 
                elif self.weir_type == 24:
                    BARINCFSB = 1.0; BARLANCFSP = 1.0
                    f_out.write('%s\t%s\n' %(int(len(list_right_weirBC)), 24))
                    for i in list_right_weirBC:
                        f_out.write('%s\t%s\t%s\t%s\t%s\n' %(int(i[0]), int(i[0]+1), i[1], BARINCFSB, BARLANCFSP))
            
            if not self.flag_weir:
                if (self.num_elems > 1):
                    f_out.write('%s\t%s\n' %(int(len(list_landleftBC)+2), 20))
                    f_out.write('%s\n' %(list_channelND[-1]))
                else:
                    f_out.write('%s\t%s\n' %(int(len(list_landleftBC)), 20))
                for i in reversed(list_landleftBC):
                    f_out.write('%s\n' %(i))
            else: 
                if self.weir_type == 23:
                    BARLANCFSP = 1.0
                    f_out.write('%s\t%s\n' %(int(len(list_left_weirBC)), 23))
                    for i in reversed(list_left_weirBC):
                        f_out.write('%s\t%s\t%s\n' %(int(i[0]), i[1], BARLANCFSP))    
                elif self.weir_type == 24:
                    BARINCFSB = 1.0; BARLANCFSP = 1.0
                    f_out.write('%s\t%s\n' %(int(len(list_left_weirBC)), 24))
                    for i in reversed(list_left_weirBC):
                        f_out.write('%s\t%s\t%s\t%s\t%s\n' %(int(i[0]), int(i[0]-1), i[1], BARINCFSB, BARLANCFSP))

        
if __name__ == '__main__':
    A = NHDFlowline()
    A.main()  
           