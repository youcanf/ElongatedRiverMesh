#// Sort and Trace Downstream
#// Author: Youcan Feng
#// Create Date: 5/25/2020
#// Purpose: Reorder the segments following the downstream direction


import os
import numpy as np
import pandas as pd
import ogr


                   
def main():

    file_path = 'Inputs/NHD_XS_pts_split2_joined.shp'
    outShapefile = "Outputs/NHD_XS_pts_split2_joined_ranked.shp" 

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(file_path,0)
    layer = dataSource.GetLayer()
    spatialRef = layer.GetSpatialRef()
    column_names = ["Fid", "StartX","StartY", "EndX", "EndY", "FromInd", "NextInd", "SortID"]
    df = pd.DataFrame(columns=column_names)

    for inFeature in layer:
        lines = inFeature.GetGeometryRef()
        fid = int(inFeature.GetFID())
        startX = lines.GetX(0)
        startY = lines.GetY(0)
        endX = lines.GetX(lines.GetPointCount()-1)
        endY = lines.GetY(lines.GetPointCount()-1)
        dict_tmp = {"Fid": fid, "StartX": startX, "StartY": startY, "EndX": endX, "EndY": endY}
        df = df.append(dict_tmp, ignore_index=True)

    for index, row in df.iterrows():
        from_tmp = df.loc[(df['EndX']==row['StartX']) & (df['EndY']==row['StartY'])]
        if len(from_tmp) == 0:
            row['FromInd'] = -9999  
        else:
            row['FromInd'] = from_tmp.index[0]
        to_tmp = df.loc[(df['StartX']==row['EndX']) & (df['StartY']==row['EndY'])]
        if len(to_tmp) == 0:
            row['NextInd'] = -9999  
        else: 
            row['NextInd'] = to_tmp.index[0]
        
    origin = df.loc[df['FromInd'] == -9999]
    if len(origin) == 1: 
        oldfid = int(origin.index[0])
    else:
        print('there is no origin.')
        exit(1)
    for i in range(len(df)):
        df.iloc[oldfid]['SortID'] = i
        oldfid = int(df.iloc[oldfid]['NextInd'])

    df['Fid'] = df['Fid'].astype(np.int64)
    df['FromInd'] = df['FromInd'].astype(np.int64)
    df['NextInd'] = df['NextInd'].astype(np.int64)
    df['SortID'] = df['SortID'].astype(np.int64)
     
    outDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outShapefile):
        outDriver.DeleteDataSource(outShapefile)
    outDataSource = outDriver.CreateDataSource(outShapefile)
    outLayer = outDataSource.CreateLayer("Rivers", spatialRef, geom_type = ogr.wkbMultiLineString)

    inLayerDefn = layer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    sortField = ogr.FieldDefn("SortID", ogr.OFTInteger64)
    outLayer.CreateField(sortField)
    outLayerDefn = outLayer.GetLayerDefn()
    
    layer.ResetReading() 
    count = 0
    for inFeature in layer:
        outFeature = ogr.Feature(outLayerDefn)
        geom = inFeature.GetGeometryRef()
        outFeature.SetGeometry(geom)
        for i in range(0,outLayerDefn.GetFieldCount()):
            if i != outLayerDefn.GetFieldCount()-1:
                outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))  
            else:
                outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), df.iloc[count]['SortID'])    
        outLayer.CreateFeature(outFeature)
        outFeature = None 
        count = count + 1            
    i = 0
  
if __name__ == '__main__':
    main()           