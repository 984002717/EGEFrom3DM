from ast import Str
from copy import deepcopy
import os
from re import T
from xmlrpc.client import Fault
import numpy as np
# from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from shapely.geometry import Polygon,MultiPoint,MultiPolygon,Point,LineString,LinearRing,MultiLineString
from shapely.ops import voronoi_diagram
import pandas as pd
import math 
if __name__=='__main__':
    FolderPath = os.path.dirname(__file__).replace('\\','/')
    StratigraphicColumn = np.loadtxt(FolderPath + "/stratigraphic_column.dat",dtype='str',delimiter='\t',skiprows=1)
    Model_Bou = np.loadtxt(FolderPath+'/bou.txt',dtype='float',delimiter='\t',skiprows=1)
    Model_Bou_line = LinearRing(Model_Bou)
    file_all = os.listdir(FolderPath)
    obj_list = []
    Tol_order = 4
    Tolerance = 10**(-Tol_order)
    Bottom = round(Model_Bou[0][2],Tol_order)
    Surface_All = {}
    Fault9nodes = []
    Fault3Id = []
    Intrusive9node = []
    Intrusive3Id = []
    fault_name = ''
    Intrusive_name = ''
    for obji in file_all:
        if os.path.splitext(obji)[1] == '.obj':
            obj_list.append(obji)
            Surface_All[obji.replace('.obj','')]={}
            Surface_All[obji.replace('.obj','')]['nodes']=[]
            Surface_All[obji.replace('.obj','')]['tins']=[]
            Surface_All[obji.replace('.obj','')]['xyz9']=[]
            Surface_All[obji.replace('.obj','')]['id3']=[]
            if obji.find('fault')>=0:
                fault_name = obji.replace('.obj','')
            if obji.lower().find('intrusive')>=0:
                Intrusive_name = obji.replace('.obj','')
    for obji in obj_list:
        f_in = open(FolderPath+'/'+obji,'r')
        node_count_read = 0
        while 1:
            lc1 = f_in.readline()
            if lc1:
                # line_id += 1
                l2 = lc1.strip().replace(" ",'\t').replace(',','\t').replace(';','\t').replace('\t\t','\t').replace('\t\t','\t').split('\t')
                if l2:
                    if l2[0].lower() == 'f':
                        if l2[1].find('/') == -1:
                            Surface_All[obji.replace('.obj','')]['tins'].append([int(l2[1]),int(l2[2]),int(l2[3])])
                        else:
                            l2[1]=l2[1][:l2[1].find('/')]
                            l2[2]=l2[2][:l2[2].find('/')]
                            l2[3]=l2[3][:l2[3].find('/')]
                            Surface_All[obji.replace('.obj','')]['tins'].append([int(l2[1]),int(l2[2]),int(l2[3])])
                    elif l2[0].lower() == 'v':
                        node_count_read += 1
                        Surface_All[obji.replace('.obj','')]['nodes'].append([round(float(l2[1]),Tol_order),round(float(l2[2]),Tol_order),round(float(l2[3]),Tol_order),int(node_count_read)])#,int(l2[1])
            else:
                break
        f_in.close()
        print(obji.replace('.obj','')+'  Read_OK!!!')
    #Get Terrian and Eroded surface of the 3D geological models
    contact_all = StratigraphicColumn[:,3]
    # count1 = contact_all
    disconformity_count = 0
    for contacti in contact_all:
        if contacti.find('disconformable')>=0 or contacti.find('unconformable')>=0:
            disconformity_count += 1
    Erod_part_count = disconformity_count+1
    Erod_part = []
    for Erodi in range(Erod_part_count):
        Erod_part.append({})
        Erod_part[Erodi]['nodes'] = []
        Erod_part[Erodi]['tins'] = []
        Erod_part[Erodi]['xyz9']=[]
        Erod_part[Erodi]['id3']=[]
    Terrian_surface = {}
    Terrian_surface['nodes'] = []
    Terrian_surface['tins'] = []
    Terrian_surface['xyz9'] = []
    Terrian_surface['id3'] = []
    Current_part = 0
    Result_Terrian = {}
    Result_Terrian['xyz9'] = []
    Result_Terrian['id3'] = []
    for layi in StratigraphicColumn:
        # obji_name = obji.replace('.obj','')
        lay_name = layi[0]
        nodes_count1 = len(Terrian_surface['nodes'])
        if lay_name+'.obj' not in obj_list:
            Surface_All[lay_name]={}
            Surface_All[lay_name]['nodes']=[]
            Surface_All[lay_name]['tins']=[]
            Surface_All[lay_name]['xyz9']=[]
            Surface_All[lay_name]['id3']=[]
        if Surface_All[lay_name]['nodes']:
            Terrian_surface['nodes'].extend(Surface_All[lay_name]['nodes'])
            Terrian_surface['tins'].extend(list(np.array(Surface_All[lay_name]['tins'])+nodes_count1))
            part_count = len(Erod_part[Current_part]['nodes'])
            Erod_part[Current_part]['nodes'].extend(Surface_All[lay_name]['nodes'])
            Erod_part[Current_part]['tins'].extend(list(np.array(Surface_All[lay_name]['tins'])+part_count))
            for tri_i in Surface_All[lay_name]['tins']:
                idd1,idd2,idd3 = tri_i
                xx1,yy1,zz1,pid1 = Surface_All[lay_name]['nodes'][idd1-1]
                xx2,yy2,zz2,pid2 = Surface_All[lay_name]['nodes'][idd2-1]
                xx3,yy3,zz3,pid3 = Surface_All[lay_name]['nodes'][idd3-1]
                dis1 = Model_Bou_line.distance(Point([xx1,yy1]))
                dis2 = Model_Bou_line.distance(Point([xx2,yy2]))
                dis3 = Model_Bou_line.distance(Point([xx3,yy3]))
                if dis1>Tolerance or dis2>Tolerance or dis3>Tolerance:
                    if abs(zz1-Bottom)>Tolerance or abs(zz2-Bottom)>Tolerance or abs(zz3-Bottom)>Tolerance:
                        delx1,dely1,delz1 = [xx2-xx1,yy2-yy1,zz2-zz1]
                        delx2,dely2,delz2 = [xx3-xx1,yy3-yy1,zz3-zz1]
                        Vetor_z = abs((delx1 * dely2) - (dely1 * delx2))
                        if Vetor_z<0.01:
                            if dis1>Tolerance*10 or dis2>Tolerance*10 or dis3>Tolerance*10:
                                Surface_All[lay_name]['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))
                                Surface_All[lay_name]['id3'].append([idd1,idd2,idd3])
                                Terrian_surface['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))#,
                                Terrian_surface['id3'].append([idd1+nodes_count1,idd2+nodes_count1,idd3+nodes_count1])
                                Erod_part[Current_part]['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))
                                Erod_part[Current_part]['id3'].append([idd1+part_count,idd2+part_count,idd3+part_count])
                        else:
                            Surface_All[lay_name]['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))
                            Surface_All[lay_name]['id3'].append([idd1,idd2,idd3])
                            Terrian_surface['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))#,
                            Terrian_surface['id3'].append([idd1+nodes_count1,idd2+nodes_count1,idd3+nodes_count1])
                            Erod_part[Current_part]['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))
                            Erod_part[Current_part]['id3'].append([idd1+part_count,idd2+part_count,idd3+part_count])
                else:
                    delx1,dely1,delz1 = [xx2-xx1,yy2-yy1,zz2-zz1]
                    delx2,dely2,delz2 = [xx3-xx1,yy3-yy1,zz3-zz1]
                    Vetor_z = abs((delx1 * dely2) - (dely1 * delx2))
                    if Vetor_z>0.01:
                    # if dis1>Tolerance*10 or dis2>Tolerance*10 or dis3>Tolerance*10:
                        Surface_All[lay_name]['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))
                        Surface_All[lay_name]['id3'].append([idd1,idd2,idd3])
                        Terrian_surface['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))#,
                        Terrian_surface['id3'].append([idd1+nodes_count1,idd2+nodes_count1,idd3+nodes_count1])
                        Erod_part[Current_part]['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))
                        Erod_part[Current_part]['id3'].append([idd1+part_count,idd2+part_count,idd3+part_count])
        if layi[3]=='disconformable' or layi[3]=='unconformable':
            Current_part += 1
    #Get Terrian Surface
    Ter_tins_count = len(Terrian_surface['xyz9'])
    nodes_count1 = len(Terrian_surface['nodes'])
    # tins_count1 = len()
    if fault_name:
        for ftri_i in Surface_All[fault_name]['tins']:
            idd1,idd2,idd3 = ftri_i
            xx1,yy1,zz1,pid1 = Surface_All[fault_name]['nodes'][idd1]
            xx2,yy2,zz2,pid2 = Surface_All[fault_name]['nodes'][idd2]
            xx3,yy3,zz3,pid3 = Surface_All[fault_name]['nodes'][idd3]
            Fault9nodes.append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))#,
            Fault3Id.append([idd1,idd2,idd3])
            # Fault3Id.append([idd1,idd2,idd3])
            Terrian_surface['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))#,
            Terrian_surface['id3'].append([idd1+nodes_count1,idd2+nodes_count1,idd3+nodes_count1])
            # f_nodes_count = len(Fault9nodes)
            # Terrian_surface['id3'].append([idd1+f_nodes_count,idd2+f_nodes_count,idd3+f_nodes_count])
    nodes_count1 = len(Terrian_surface['nodes'])
    if Intrusive_name:
        for ftri_i in Surface_All[Intrusive_name]['tins']:
            idd1,idd2,idd3 = ftri_i
            xx1,yy1,zz1,pid1 = Surface_All[Intrusive_name]['nodes'][idd1-1]
            xx2,yy2,zz2,pid2 = Surface_All[Intrusive_name]['nodes'][idd2-1]
            xx3,yy3,zz3,pid3 = Surface_All[Intrusive_name]['nodes'][idd3-1]
            Intrusive9node.append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))#,
            Intrusive3Id.append([idd1,idd2,idd3])
            Terrian_surface['xyz9'].append(sorted([xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3]))#,
            Terrian_surface['id3'].append([idd1+nodes_count1,idd2+nodes_count1,idd3+nodes_count1])
            # i_nodes_count = len(Fault9nodes)
            # Terrian_surface['id3'].append([idd1+i_nodes_count,idd2+i_nodes_count,idd3+i_nodes_count])

    NodesAll_faultIntru = Terrian_surface['xyz9']
    NodesIdAll_faultTintru = Terrian_surface['id3']
    # if Fault9nodes:
    #     NodesAll_faultIntru.extend(Fault9nodes)
    #     NodesIdAll_faultTintru.extend(Fault3Id)
    # if Intrusive9node:
    #     NodesAll_faultIntru.extend(Intrusive9node)
    #     NodesIdAll_faultTintru.append(Intrusive3Id)
    Ter_tinsNodes_np = np.array(NodesAll_faultIntru)
    
    # Ter_tinsId_np = np.array(Terrian_surface['id3'])
    Ter_tinsId_np = np.array(NodesIdAll_faultTintru)
    ter_pds = pd.DataFrame({
        'xx1':list(Ter_tinsNodes_np[:,0].flatten()),
        'yy1':list(Ter_tinsNodes_np[:,1].flatten()),
        'zz1':list(Ter_tinsNodes_np[:,2].flatten()),
        'xx2':list(Ter_tinsNodes_np[:,3].flatten()),
        'yy2':list(Ter_tinsNodes_np[:,4].flatten()),
        'zz2':list(Ter_tinsNodes_np[:,5].flatten()),
        'xx3':list(Ter_tinsNodes_np[:,6].flatten()),
        'yy3':list(Ter_tinsNodes_np[:,7].flatten()),
        'zz3':list(Ter_tinsNodes_np[:,8].flatten())
    })
    same1 = ter_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'last')#保留最后一个，，其他的标记True
    same2 = ter_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'first')#除第一个其他标记True
    same_all = np.array(list(zip(list(same1),list(same2))))
    #same_all = np.hstack((same_all[np.arange(Ter_tins_count)],np.ones((same_all.shape[0]-Ter_tins_countc,2),dtype=bool)))
    same_all_result = same_all.any(axis = 1)#获取所有的True，后面提取这些重复的点
    trgl_same_reverse = np.logical_not(same_all_result)
    trgl_same_reverse = np.hstack((trgl_same_reverse[np.arange(Ter_tins_count)],np.zeros(trgl_same_reverse.shape[0]-Ter_tins_count,dtype=bool)))
    trgl_result = Ter_tinsId_np[trgl_same_reverse]
    Result_Terrian['xyz9'] = list(Ter_tinsNodes_np[trgl_same_reverse])
    Result_Terrian['id3'] = list(Ter_tinsId_np[trgl_same_reverse])

    all_related_nodesId = np.unique(np.sort(trgl_result.flatten()))
    ter_id = {}
    for teri in range(all_related_nodesId.shape[0]):
        ter_id[str(all_related_nodesId[teri])] = teri+1
    all_related_nodes = np.array(Terrian_surface['nodes'])[all_related_nodesId-1]
    Result_Terrian['nodes'] = deepcopy(all_related_nodes)
    if not os.path.exists(FolderPath+'/Terrian'):
        os.mkdir(FolderPath+'/Terrian')
    f_out = open(FolderPath+'/Terrian/Dem.obj','w')
    f_out.write('g TerrianResult\n')
    for nid in all_related_nodes:
        if abs(nid[0]--842.552001953125)<0.00001 and abs(nid[1]-97.998054504394531)<0.0001 and abs(nid[2]-60.899604797363281)<0.0001:
            fggrggrgrrg=1
        f_out.write('v '+ str(nid[0])+' '+str(nid[1])+' '+str(nid[2])+'\n')
    for nid in trgl_result:
        f_out.write('f '+ str(ter_id[str(nid[0])])+' '+str(ter_id[str(nid[1])])+' '+str(ter_id[str(nid[2])])+'\n')
    print('Get Terrian surface OK...')
    # Get n Eroded surface of the model
    Result_Eroded_Surface = []
    for ero_i in range(Erod_part_count-1,0,-1):
        if ero_i == 0:
            continue
        else:
            Result_Eroded_Surface.append({})
            Result_Eroded_Surface[-1]['nodes']=[]
            Result_Eroded_Surface[-1]['tins']=[]
            Result_Eroded_Surface[-1]['xyz9']=[]
            Result_Eroded_Surface[-1]['id3']=[]
            Ero_i_tinNodes_np = Erod_part[ero_i]['xyz9']
            Ero_i_tinId_np = Erod_part[ero_i]['id3']
            Erod_i_Tins_count = len(Erod_part[ero_i]['xyz9'])
            if Erod_part[ero_i]['xyz9']:
                if fault_name:#Add Fault Tins to current Eroded surface
                    Ero_i_tinNodes_np.extend(Fault9nodes)
                    Ero_i_tinId_np.extend(Fault3Id)
                if Intrusive_name:#Add IntrusiveBody Tins to current Eroded surface
                    Ero_i_tinNodes_np.extend(Intrusive9node)
                    Ero_i_tinId_np.extend(Intrusive3Id)
                for beforeEroi in range(0,len(Result_Eroded_Surface)):#Add Eroded suface, formed before current Eroded, Tins to current Eroded surface
                    Ero_i_tinNodes_np.extend(Result_Eroded_Surface[beforeEroi]['xyz9'])
                    Ero_i_tinId_np.extend(Result_Eroded_Surface[beforeEroi]['id3'])
                Ero_i_tinNodes_np = np.array(Ero_i_tinNodes_np)
                Ero_i_tinId_np = np.array(Ero_i_tinId_np)
                # Result_Eroded_Surface[-1]['id3']
                Ero_i_pds = pd.DataFrame({
                    'xx1':list(Ero_i_tinNodes_np[:,0].flatten()),
                    'yy1':list(Ero_i_tinNodes_np[:,1].flatten()),
                    'zz1':list(Ero_i_tinNodes_np[:,2].flatten()),
                    'xx2':list(Ero_i_tinNodes_np[:,3].flatten()),
                    'yy2':list(Ero_i_tinNodes_np[:,4].flatten()),
                    'zz2':list(Ero_i_tinNodes_np[:,5].flatten()),
                    'xx3':list(Ero_i_tinNodes_np[:,6].flatten()),
                    'yy3':list(Ero_i_tinNodes_np[:,7].flatten()),
                    'zz3':list(Ero_i_tinNodes_np[:,8].flatten())
                })
                same11 = Ero_i_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'last')#保留最后一个，，其他的标记True
                same22 = Ero_i_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'first')#除第一个其他标记True
                same_all2 = np.array(list(zip(list(same11),list(same22))))
                same_all_result2 = same_all2.any(axis = 1)#获取所有的True，后面提取这些重复的点
                trgl_same_reverse2 = np.logical_not(same_all_result2)
                trgl_same_reverse2 = np.hstack((trgl_same_reverse2[np.arange(Erod_i_Tins_count)],np.zeros(trgl_same_reverse2.shape[0]-Erod_i_Tins_count,dtype=bool)))
                trgl_result2 = Ero_i_tinId_np[trgl_same_reverse2]
                Result_Eroded_Surface[-1]['xyz9'] = list(Ero_i_tinNodes_np[trgl_same_reverse2])
                Result_Eroded_Surface[-1]['id3'] = list(Ero_i_tinId_np[trgl_same_reverse2])
                Ero_i_all_related_nodesId = np.unique(np.sort(trgl_result2.flatten()))
                ter_id = {}
                for teri in range(Ero_i_all_related_nodesId.shape[0]):
                    ter_id[str(Ero_i_all_related_nodesId[teri])] = teri+1
                Ero_i_all_related_nodes = np.array(Erod_part[ero_i]['nodes'])[Ero_i_all_related_nodesId-1]
                Result_Eroded_Surface[-1]['nodes'] = deepcopy(list(Ero_i_all_related_nodes))
                if not os.path.exists(FolderPath+'/ErodedSurface'):
                    os.mkdir(FolderPath+'/ErodedSurface')
                f_out = open(FolderPath+'/ErodedSurface/ErodedSurface_'+str(len(Result_Eroded_Surface))+'.obj','w')
                f_out.write('g ErodedSurface_'+str(len(Result_Eroded_Surface))+'\n')
                for nid in Ero_i_all_related_nodes:
                    Result_Eroded_Surface[-1]['nodes'].append([nid[0],nid[1],nid[2],len(Result_Eroded_Surface[-1]['nodes'])+1])
                    f_out.write('v '+ str(nid[0])+' '+str(nid[1])+' '+str(nid[2])+'\n')
                for nid in trgl_result2:
                    f_out.write('f '+ str(ter_id[str(nid[0])])+' '+str(ter_id[str(nid[1])])+' '+str(ter_id[str(nid[2])])+'\n')
                    Result_Eroded_Surface[-1]['tins'].append([ter_id[str(nid[0])],ter_id[str(nid[1])],ter_id[str(nid[2])]])
                f_out.close()
    print('Get '+str(Erod_part_count-1)+' Eroded surface OK...')
    #Get Top and Bottom Surface
    Result_top_Surface = {}
    Result_bot_Surface = {}
    Result_top_bou_Surface = {}
    Result_bot_bou_Surface = {}
    for layi in range(StratigraphicColumn.shape[0]):
        # obji_name = obji.replace('.obj','')
        lay_name = StratigraphicColumn[layi][0]
        Result_top_Surface[lay_name]={}
        Result_top_Surface[lay_name]['nodes']=[]
        Result_top_Surface[lay_name]['tins']=[]
        Result_top_Surface[lay_name]['xyz9']=[]
        Result_top_Surface[lay_name]['id3']=[]
        Result_bot_Surface[lay_name] = {}
        Result_bot_Surface[lay_name]['nodes']=[]
        Result_bot_Surface[lay_name]['tins']=[]
        Result_bot_Surface[lay_name]['xyz9']=[]
        Result_bot_Surface[lay_name]['id3']=[]
        # For Top Suface
        cur_sur_tinsNode_count = len(Surface_All[lay_name]['xyz9'])
        cur_sur_top_tinsNode_np = deepcopy(Surface_All[lay_name]['xyz9'])
        cur_sur_top_tinsId_np = deepcopy(Surface_All[lay_name]['id3'])
        cur_sur_bot_tinsNode_np = deepcopy(Surface_All[lay_name]['xyz9'])
        cur_sur_bot_tinsId_np = deepcopy(Surface_All[lay_name]['id3'])
        if fault_name:
            cur_sur_top_tinsNode_np.extend(Fault9nodes)
            cur_sur_top_tinsId_np.extend(Fault3Id)
            cur_sur_bot_tinsNode_np.extend(Fault9nodes)
            cur_sur_bot_tinsId_np.extend(Fault3Id)
        if Intrusive_name:
            cur_sur_top_tinsNode_np.extend(Intrusive9node)
            cur_sur_top_tinsId_np.extend(Intrusive3Id)
            cur_sur_bot_tinsNode_np.extend(Intrusive9node)
            cur_sur_bot_tinsId_np.extend(Intrusive3Id)
        if Erod_part_count > 1:
            for beforeEroi in range(0,len(Result_Eroded_Surface)):#Add Eroded suface, formed before current Eroded, Tins to current Eroded surface
                cur_sur_top_tinsNode_np.extend(Result_Eroded_Surface[beforeEroi]['xyz9'])
                cur_sur_top_tinsId_np.extend(Result_Eroded_Surface[beforeEroi]['id3'])
        after_Eroded_count = 0
        for lli in range(layi):
            if StratigraphicColumn[lli][3] == 'disconformable' or StratigraphicColumn[lli][3]=='unconformable':
                after_Eroded_count += 1
        for erii in range(after_Eroded_count):
            cur_sur_bot_tinsNode_np.extend(Result_Eroded_Surface[erii]['xyz9'])
            cur_sur_bot_tinsId_np.extend(Result_Eroded_Surface[erii]['id3'])
        # lay_below_cur = StratigraphicColumn[index_id+1:,:]
        for lay_low_i in range(layi+1,StratigraphicColumn.shape[0]):
            below_name = StratigraphicColumn[lay_low_i][0]
            cur_sur_top_tinsNode_np.extend(Surface_All[below_name]['xyz9'])
            cur_sur_top_tinsId_np.extend(Surface_All[below_name]['id3'])
        #if lay_up_cur:
        for lay_up_i in range(layi):
            up_name = StratigraphicColumn[lay_up_i][0]
            cur_sur_bot_tinsNode_np.extend(Surface_All[up_name]['xyz9'])
            cur_sur_bot_tinsId_np.extend(Surface_All[up_name]['id3'])
        # #Add Terrian Surface to current surface
        # cur_sur_bot_tinsNode_np.extend(Result_Terrian['xyz9'])
        # cur_sur_top_tinsNode_np.extend(Result_Terrian['xyz9'])
        # cur_sur_bot_tinsId_np.extend(Result_Terrian['id3'])
        # cur_sur_top_tinsId_np.extend(Result_Terrian['id3'])
        cur_sur_top_tinsNode_np = np.array(cur_sur_top_tinsNode_np)
        cur_sur_top_tinsId_np = np.array(cur_sur_top_tinsId_np)
        cur_sur_pds = pd.DataFrame({
            'xx1':list(cur_sur_top_tinsNode_np[:,0].flatten()),
            'yy1':list(cur_sur_top_tinsNode_np[:,1].flatten()),
            'zz1':list(cur_sur_top_tinsNode_np[:,2].flatten()),
            'xx2':list(cur_sur_top_tinsNode_np[:,3].flatten()),
            'yy2':list(cur_sur_top_tinsNode_np[:,4].flatten()),
            'zz2':list(cur_sur_top_tinsNode_np[:,5].flatten()),
            'xx3':list(cur_sur_top_tinsNode_np[:,6].flatten()),
            'yy3':list(cur_sur_top_tinsNode_np[:,7].flatten()),
            'zz3':list(cur_sur_top_tinsNode_np[:,8].flatten())
        })
        same11 = cur_sur_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'last')#保留最后一个，，其他的标记True
        same22 = cur_sur_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'first')#除第一个其他标记True
        same_all2 = np.array(list(zip(list(same11),list(same22))))
        same_all_result2 = same_all2.any(axis = 1)#获取所有的True，后面提取这些重复的点
        trgl_same_reverse2 = np.logical_not(same_all_result2)
        trgl_same_reverse2 = np.hstack((trgl_same_reverse2[np.arange(cur_sur_tinsNode_count)],np.zeros(trgl_same_reverse2.shape[0]-cur_sur_tinsNode_count,dtype=bool)))
        trgl_result2 = cur_sur_top_tinsId_np[trgl_same_reverse2]
        Result_top_Surface[lay_name]['xyz9'] = list(cur_sur_top_tinsNode_np[trgl_same_reverse2])
        Result_top_Surface[lay_name]['id3'] = list(cur_sur_top_tinsId_np[trgl_same_reverse2])
        Result_top_all_related_nodesId = np.unique(np.sort(trgl_result2.flatten()))
        ter_id = {}
        for teri in range(Result_top_all_related_nodesId.shape[0]):
            ter_id[str(Result_top_all_related_nodesId[teri])] = teri+1
        Result_top_all_related_nodes = np.array(Surface_All[lay_name]['nodes'])[Result_top_all_related_nodesId-1]
        if not os.path.exists(FolderPath+'/TopBotSurface'):
            os.mkdir(FolderPath+'/TopBotSurface')
        f_out = open(FolderPath+'/TopBotSurface/'+lay_name+'_Top.obj','w')
        f_out.write('g Top_'+lay_name+'\n')
        for nid in Result_top_all_related_nodes:
            Result_top_Surface[lay_name]['nodes'].append([nid[0],nid[1],nid[2],len(Result_top_Surface[lay_name]['nodes'])+1])
            f_out.write('v '+ str(nid[0])+' '+str(nid[1])+' '+str(nid[2])+'\n')
        for nid in trgl_result2:
            f_out.write('f '+ str(ter_id[str(nid[0])])+' '+str(ter_id[str(nid[1])])+' '+str(ter_id[str(nid[2])])+'\n')
            Result_top_Surface[lay_name]['tins'].append([ter_id[str(nid[0])],ter_id[str(nid[1])],ter_id[str(nid[2])]])
        f_out.close()
        # Extract Top surface Boundary
        if trgl_result2.size>0:
            TopBou = MultiLineString()
            trgl_result2_p1 = trgl_result2[:,0]
            trgl_result2_p2 = trgl_result2[:,1]
            trgl_result2_p3 = trgl_result2[:,2]
            nodes_col_1 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p1-1][:,:3]
            nodes_col_2 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p2-1][:,:3]
            nodes_col_3 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p3-1][:,:3]
            Line1 = np.hstack((nodes_col_1,nodes_col_2))
            Line2 = np.hstack((nodes_col_1,nodes_col_3))
            Line3 = np.hstack((nodes_col_2,nodes_col_3))
            LineAll = np.vstack((Line1,Line2))
            LineAll = np.vstack((LineAll,Line3))
            Line_init = deepcopy(LineAll)
            LineAll.sort(axis=1)
            Top_Bou_pd = pd.DataFrame({
                'xx1':list(LineAll[:,0].flatten()),
                'yy1':list(LineAll[:,1].flatten()),
                'zz1':list(LineAll[:,2].flatten()),
                'xx2':list(LineAll[:,3].flatten()),
                'yy2':list(LineAll[:,4].flatten()),
                'zz2':list(LineAll[:,5].flatten()),
            })
            same11 = Top_Bou_pd.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2'],keep = 'last')#保留最后一个，，其他的标记True
            same22 = Top_Bou_pd.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2'],keep = 'first')#除第一个其他标记True
            same_all2 = np.array(list(zip(list(same11),list(same22))))
            same_all_result2 = same_all2.any(axis = 1)#获取所有的True，后面提取这些重复的点
            trgl_same_reverse2 = np.logical_not(same_all_result2)
            Line_part = Line_init[trgl_same_reverse2]
            if Line_part.size>0:
                if Line_part.ndim ==2:
                    for l_i in Line_part:
                        xx1,yy1,zz1,xx2,yy2,zz2 = l_i
                        TopBou = TopBou.union(LineString([[xx1,yy1,zz1],[xx2,yy2,zz2]]))
            if TopBou:
                f_out2 = open(FolderPath+'/TopBotSurface/'+lay_name+'_TopBou.dat','w')
                f_out3 = open(FolderPath+'/TopBotSurface/'+lay_name+'_TopBou_ContactTer.dat','w')
                f_out4 = open(FolderPath+'/TopBotSurface/'+lay_name+'_TopBou_ContactErod.dat','w')
                Terrian_Nodes_all = np.array(Result_Terrian['nodes'])[:,:3].tolist()
                Erod_Nodes_all = []
                for erri in range(len(Result_Eroded_Surface)):
                    Erod_Nodes_all.extend(np.array(Result_Eroded_Surface[erri]['nodes'])[:,:3].tolist())
                if TopBou.geom_type == "LineString" or TopBou.geom_type == "LinearRing":
                    f_out2.write('line 1 coords\n')
                    coords_bou = np.array(TopBou.wkt[TopBou.wkt.find('(')+1:-1].replace(', ',',').replace(' ',',').split(','))
                    coords_bou = coords_bou.reshape((-1,3)).astype('float')
                    coord_bou_p_in_Ter = []
                    coord_bou_p_in_Erod = []
                    for bou_coord in coords_bou:
                        f_out2.write(str(bou_coord[0])+','+str(bou_coord[1])+','+str(bou_coord[2])+'\n')
                        if [bou_coord[0],bou_coord[1],bou_coord[2]] in Terrian_Nodes_all:
                            coord_bou_p_in_Ter.append(True)
                        else:
                            coord_bou_p_in_Ter.append(False)
                        if [bou_coord[0],bou_coord[1],bou_coord[2]] in Erod_Nodes_all:
                            coord_bou_p_in_Erod.append(True)
                        else:
                            coord_bou_p_in_Erod.append(False)
                    bou_coord_count = coords_bou.shape[0]
                    first_in_Ter = coord_bou_p_in_Ter[0]
                    if first_in_Ter == True:
                        f_out3.write('line\n')
                    first_in_Erod = coord_bou_p_in_Erod[0]
                    if first_in_Erod == True:
                        f_out4.write('line\n')
                    for bou_coord in range(bou_coord_count):
                        nextii = (bou_coord+1) % bou_coord_count
                        preii = (bou_coord-1+bou_coord_count)%bou_coord_count
                        if coord_bou_p_in_Ter[bou_coord] == True:
                            if coord_bou_p_in_Ter[nextii] == True:
                                if coord_bou_p_in_Ter[preii] ==False:
                                    f_out3.write('line\n')
                                f_out3.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                            else:
                                if coord_bou_p_in_Ter[preii] == True:
                                    f_out3.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                        if coord_bou_p_in_Erod[bou_coord] == True:
                            if coord_bou_p_in_Erod[nextii] == True:
                                if coord_bou_p_in_Erod[preii] ==False:
                                    f_out4.write('line\n')
                                f_out4.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                            else:
                                if coord_bou_p_in_Erod[preii] ==True:
                                    f_out4.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                
                elif TopBou.geom_type == "MultiLineString":
                    for boui in range(len(TopBou.geoms)):
                        bou_temp = TopBou.geoms[boui]
                        f_out2.write('line '+str(boui)+' coords\n')
                        coords_bou = np.array(bou_temp.wkt[bou_temp.wkt.find('(')+1:-1].replace(', ',',').replace(' ',',').split(','))
                        coords_bou = coords_bou.reshape((-1,3)).astype('float')
                        coord_bou_p_in_Ter = []
                        coord_bou_p_in_Erod = []
                        for bou_coord in coords_bou:
                            f_out2.write(str(bou_coord[0])+','+str(bou_coord[1])+','+str(bou_coord[2])+'\n')
                            if [bou_coord[0],bou_coord[1],bou_coord[2]] in Terrian_Nodes_all:
                                coord_bou_p_in_Ter.append(True)
                            else:
                                coord_bou_p_in_Ter.append(False)
                            if [bou_coord[0],bou_coord[1],bou_coord[2]] in Erod_Nodes_all:
                                coord_bou_p_in_Erod.append(True)
                            else:
                                coord_bou_p_in_Erod.append(False)
                        bou_coord_count = coords_bou.shape[0]
                        first_in_Ter = coord_bou_p_in_Ter[0]
                        if first_in_Ter == True:
                            f_out3.write('line\n')
                        first_in_Erod = coord_bou_p_in_Erod[0]
                        if first_in_Erod == True:
                            f_out4.write('line\n')
                        for bou_coord in range(bou_coord_count):
                            nextii = (bou_coord+1) % bou_coord_count
                            preii = (bou_coord-1+bou_coord_count)%bou_coord_count
                            if coord_bou_p_in_Ter[bou_coord] == True:
                                if coord_bou_p_in_Ter[nextii] == True:
                                    if coord_bou_p_in_Ter[preii] ==False:
                                        f_out3.write('line\n')
                                    f_out3.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                                else:
                                    if coord_bou_p_in_Ter[preii] == True:
                                        f_out3.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                            if coord_bou_p_in_Erod[bou_coord] == True:
                                if coord_bou_p_in_Erod[nextii] == True:
                                    if coord_bou_p_in_Erod[preii] ==False:
                                        f_out4.write('line\n')
                                    f_out4.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                                else:
                                    if coord_bou_p_in_Erod[preii] ==True:
                                        f_out4.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                f_out3.close()
                f_out4.close()
                f_out2.close()
        #Add Top Surafce to Current Bottom surface
        cur_sur_bot_tinsNode_np.extend(Result_top_Surface[lay_name]['xyz9'])# = list(cur_sur_top_tinsNode_np[trgl_same_reverse2])
        cur_sur_bot_tinsId_np.extend(Result_top_Surface[lay_name]['id3'])# = list(cur_sur_top_tinsId_np[trgl_same_reverse2])
        cur_sur_bot_tinsNode_np = np.array(cur_sur_bot_tinsNode_np)
        cur_sur_bot_tinsId_np = np.array(cur_sur_bot_tinsId_np)
        cur_sur_bot_pds = pd.DataFrame({
            'xx1':list(cur_sur_bot_tinsNode_np[:,0].flatten()),
            'yy1':list(cur_sur_bot_tinsNode_np[:,1].flatten()),
            'zz1':list(cur_sur_bot_tinsNode_np[:,2].flatten()),
            'xx2':list(cur_sur_bot_tinsNode_np[:,3].flatten()),
            'yy2':list(cur_sur_bot_tinsNode_np[:,4].flatten()),
            'zz2':list(cur_sur_bot_tinsNode_np[:,5].flatten()),
            'xx3':list(cur_sur_bot_tinsNode_np[:,6].flatten()),
            'yy3':list(cur_sur_bot_tinsNode_np[:,7].flatten()),
            'zz3':list(cur_sur_bot_tinsNode_np[:,8].flatten())
        })
        same11 = cur_sur_bot_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'last')#保留最后一个，，其他的标记True
        same22 = cur_sur_bot_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'first')#除第一个其他标记True
        same_all2 = np.array(list(zip(list(same11),list(same22))))
        same_all_result2 = same_all2.any(axis = 1)#获取所有的True，后面提取这些重复的点
        trgl_same_reverse2 = np.logical_not(same_all_result2)
        trgl_same_reverse2 = np.hstack((trgl_same_reverse2[np.arange(cur_sur_tinsNode_count)],np.zeros(trgl_same_reverse2.shape[0]-cur_sur_tinsNode_count,dtype=bool)))
        trgl_result2 = cur_sur_bot_tinsId_np[trgl_same_reverse2]
        Result_bot_Surface[lay_name]['xyz9'] = list(cur_sur_bot_tinsNode_np[trgl_same_reverse2])
        Result_bot_Surface[lay_name]['id3'] = list(cur_sur_bot_tinsId_np[trgl_same_reverse2])
        Result_bot_all_related_nodesId = np.unique(np.sort(trgl_result2.flatten()))
        ter_id = {}
        for teri in range(Result_bot_all_related_nodesId.shape[0]):
            ter_id[str(Result_bot_all_related_nodesId[teri])] = teri+1
        Result_bot_all_related_nodes = np.array(Surface_All[lay_name]['nodes'])[Result_bot_all_related_nodesId-1]
        if not os.path.exists(FolderPath+'/TopBotSurface'):
            os.mkdir(FolderPath+'/TopBotSurface')
        f_out = open(FolderPath+'/TopBotSurface/'+lay_name+'_Bottom.obj','w')
        f_out.write('g Bot_'+lay_name+'\n')
        for nid in Result_bot_all_related_nodes:
            Result_bot_Surface[lay_name]['nodes'].append([nid[0],nid[1],nid[2],len(Result_bot_Surface[lay_name]['nodes'])+1])
            f_out.write('v '+ str(nid[0])+' '+str(nid[1])+' '+str(nid[2])+'\n')
        for nid in trgl_result2:
            f_out.write('f '+ str(ter_id[str(nid[0])])+' '+str(ter_id[str(nid[1])])+' '+str(ter_id[str(nid[2])])+'\n')
            Result_bot_Surface[lay_name]['tins'].append([ter_id[str(nid[0])],ter_id[str(nid[1])],ter_id[str(nid[2])]])
        f_out.close()
        # Extract Bottom surface Boundary
        if trgl_result2.size>0:
            BotBou = MultiPolygon()
            trgl_result2_p1 = trgl_result2[:,0]
            trgl_result2_p2 = trgl_result2[:,1]
            trgl_result2_p3 = trgl_result2[:,2]
            nodes_col_1 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p1-1][:,:3]
            nodes_col_2 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p2-1][:,:3]
            nodes_col_3 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p3-1][:,:3]
            Line1 = np.hstack((nodes_col_1,nodes_col_2))
            Line2 = np.hstack((nodes_col_1,nodes_col_3))
            Line3 = np.hstack((nodes_col_2,nodes_col_3))
            LineAll = np.vstack((Line1,Line2))
            LineAll = np.vstack((LineAll,Line3))
            Line_init = deepcopy(LineAll)
            LineAll.sort(axis=1)
            Bot_Bou_pd = pd.DataFrame({
                'xx1':list(LineAll[:,0].flatten()),
                'yy1':list(LineAll[:,1].flatten()),
                'zz1':list(LineAll[:,2].flatten()),
                'xx2':list(LineAll[:,3].flatten()),
                'yy2':list(LineAll[:,4].flatten()),
                'zz2':list(LineAll[:,5].flatten()),
            })
            same11 = Bot_Bou_pd.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2'],keep = 'last')#保留最后一个，，其他的标记True
            same22 = Bot_Bou_pd.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2'],keep = 'first')#除第一个其他标记True
            same_all2 = np.array(list(zip(list(same11),list(same22))))
            same_all_result2 = same_all2.any(axis = 1)#获取所有的True，后面提取这些重复的点
            trgl_same_reverse2 = np.logical_not(same_all_result2)
            Line_part = Line_init[trgl_same_reverse2]
            if Line_part.size>0:
                if Line_part.ndim ==2:
                    for l_i in Line_part:
                        xx1,yy1,zz1,xx2,yy2,zz2 = l_i
                        BotBou = BotBou.union(LineString([[xx1,yy1,zz1],[xx2,yy2,zz2]]))
            if BotBou:
                f_out2 = open(FolderPath+'/TopBotSurface/'+lay_name+'_BotBou.dat','w')
                f_out3 = open(FolderPath+'/TopBotSurface/'+lay_name+'_BotBou_ContactTer.dat','w')
                f_out4 = open(FolderPath+'/TopBotSurface/'+lay_name+'_BotBou_ContactErod.dat','w')
                Terrian_Nodes_all = np.array(Result_Terrian['nodes'])[:,:3].tolist()
                Erod_Nodes_all = []
                for erri in range(len(Result_Eroded_Surface)):
                    Erod_Nodes_all.extend(np.array(Result_Eroded_Surface[erri]['nodes'])[:,:3].tolist())
                if BotBou.geom_type == "LineString" or BotBou.geom_type == "LinearRing":
                    f_out2.write('line 1 coords\n')
                    coords_bou = np.array(BotBou.wkt[BotBou.wkt.find('(')+1:-1].replace(', ',',').replace(' ',',').split(','))
                    coords_bou = coords_bou.reshape((-1,3)).astype('float')
                    coord_bou_p_in_Ter = []
                    coord_bou_p_in_Erod = []
                    for bou_coord in coords_bou:
                        f_out2.write(str(bou_coord[0])+','+str(bou_coord[1])+','+str(bou_coord[2])+'\n')
                        if [bou_coord[0],bou_coord[1],bou_coord[2]] in Terrian_Nodes_all:
                            coord_bou_p_in_Ter.append(True)
                        else:
                            coord_bou_p_in_Ter.append(False)
                        if [bou_coord[0],bou_coord[1],bou_coord[2]] in Erod_Nodes_all:
                            coord_bou_p_in_Erod.append(True)
                        else:
                            coord_bou_p_in_Erod.append(False)
                    bou_coord_count = coords_bou.shape[0]
                    first_in_Ter = coord_bou_p_in_Ter[0]
                    if first_in_Ter == True:
                        f_out3.write('line\n')
                    first_in_Erod = coord_bou_p_in_Erod[0]
                    if first_in_Erod == True:
                        f_out4.write('line\n')
                    for bou_coord in range(bou_coord_count):
                        nextii = (bou_coord+1) % bou_coord_count
                        preii = (bou_coord-1+bou_coord_count)%bou_coord_count
                        if coord_bou_p_in_Ter[bou_coord] == True:
                            if coord_bou_p_in_Ter[nextii] == True:
                                if coord_bou_p_in_Ter[preii] ==False:
                                    f_out3.write('line\n')
                                f_out3.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                            else:
                                if coord_bou_p_in_Ter[preii] == True:
                                    f_out3.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                        if coord_bou_p_in_Erod[bou_coord] == True:
                            if coord_bou_p_in_Erod[nextii] == True:
                                if coord_bou_p_in_Erod[preii] ==False:
                                    f_out4.write('line\n')
                                f_out4.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                            else:
                                if coord_bou_p_in_Erod[preii] ==True:
                                    f_out4.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                    # bou_coord_count = coords_bou.shape[0]
                    # first_in_Erod = coord_bou_p_in_Erod[0]
                    # if first_in_Erod == True:
                    #     f_out4.write('line\n')
                    # for bou_coord in range(bou_coord_count):,
                    #     nextii = (bou_coord+1) % bou_coord_count
                    #     preii = (bou_coord-1+bou_coord_count)%bou_coord_count
                    #     if coord_bou_p_in_Erod[nextii] == True and coord_bou_p_in_Erod[bou_coord] == True::
                    #         if coord_bou_p_in_Erod[preii] ==False:
                    #             f_out4.write('line\n')
                    #         f_out4.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                elif BotBou.geom_type == "MultiLineString":
                    for boui in range(len(BotBou.geoms)):
                        bou_temp = BotBou.geoms[boui]
                        f_out2.write('line '+str(boui)+' coords\n')
                        coords_bou = np.array(bou_temp.wkt[bou_temp.wkt.find('(')+1:-1].replace(', ',',').replace(' ',',').split(','))
                        coords_bou = coords_bou.reshape((-1,3)).astype('float')
                        coord_bou_p_in_Ter = []
                        coord_bou_p_in_Erod = []
                        for bou_coord in coords_bou:
                            f_out2.write(str(bou_coord[0])+','+str(bou_coord[1])+','+str(bou_coord[2])+'\n')
                            if [bou_coord[0],bou_coord[1],bou_coord[2]] in Terrian_Nodes_all:
                                coord_bou_p_in_Ter.append(True)
                            else:
                                coord_bou_p_in_Ter.append(False)
                            if [bou_coord[0],bou_coord[1],bou_coord[2]] in Erod_Nodes_all:
                                coord_bou_p_in_Erod.append(True)
                            else:
                                coord_bou_p_in_Erod.append(False)
                        bou_coord_count = coords_bou.shape[0]
                        first_in_Ter = coord_bou_p_in_Ter[0]
                        if first_in_Ter == True:
                            f_out3.write('line\n')
                        first_in_Erod = coord_bou_p_in_Erod[0]
                        if first_in_Erod == True:
                            f_out4.write('line\n')
                        for bou_coord in range(bou_coord_count):
                            nextii = (bou_coord+1) % bou_coord_count
                            preii = (bou_coord-1+bou_coord_count)%bou_coord_count
                            if coord_bou_p_in_Ter[bou_coord] == True:
                                if coord_bou_p_in_Ter[nextii] == True:
                                    if coord_bou_p_in_Ter[preii] ==False:
                                        f_out3.write('line\n')
                                    f_out3.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                                else:
                                    if coord_bou_p_in_Ter[preii] == True:
                                        f_out3.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                            if coord_bou_p_in_Erod[bou_coord] == True:
                                if coord_bou_p_in_Erod[nextii] == True:
                                    if coord_bou_p_in_Erod[preii] ==False:
                                        f_out4.write('line\n')
                                    f_out4.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                                else:
                                    if coord_bou_p_in_Erod[preii] ==True:
                                        f_out4.write(str(coords_bou[bou_coord][0])+','+str(coords_bou[bou_coord][1])+','+str(coords_bou[bou_coord][2])+'\n')
                f_out3.close()
                f_out4.close()
                f_out2.close()


        print(lay_name+' Top & Bot Surface & Bou extract OK...')
        #Get outcrop Boundary of Top and Bottom surface
        Result_top_bou_Surface[lay_name]={}
        Result_top_bou_Surface[lay_name]['nodes']=[]
        Result_top_bou_Surface[lay_name]['tins']=[]
        Result_top_bou_Surface[lay_name]['xyz9']=[]
        Result_top_bou_Surface[lay_name]['id3']=[]
        # Result_bot_bou_Surface[lay_name] = {}
        # Result_bot_bou_Surface[lay_name]['nodes']=[]
        # Result_bot_bou_Surface[lay_name]['tins']=[]
        # Result_bot_bou_Surface[lay_name]['xyz9']=[]
        # Result_bot_bou_Surface[lay_name]['id3']=[]
        cur_sur_top_bou_tinsNode_np = deepcopy(Surface_All[lay_name]['xyz9'])
        cur_sur_top_bou_tinsId_np = deepcopy(Surface_All[lay_name]['id3'])
        # cur_sur_bot_bou_tinsNode_np = deepcopy(Surface_All[lay_name]['xyz9'])
        # cur_sur_bot_bou_tinsId_np = deepcopy(Surface_All[lay_name]['id3'])
        # cur_sur_bot_tinsNode_np.extend(Result_Terrian['xyz9'])
        # cur_sur_top_tinsNode_np.extend(Result_Terrian['xyz9'])
        # cur_sur_bot_tinsId_np.extend(Result_Terrian['id3'])
        # cur_sur_top_tinsId_np.extend(Result_Terrian['id3'])
        # cur_sur_bot_tinsNode_np.extend(Result_Terrian['xyz9'])
        cur_sur_top_bou_tinsNode_np.extend(Result_Terrian['xyz9'])
        # cur_sur_bot_tinsId_np.extend(Result_Terrian['id3'])
        cur_sur_top_bou_tinsId_np.extend(Result_Terrian['id3'])
        cur_sur_top_bou_tinsNode_np = np.array(cur_sur_top_bou_tinsNode_np)
        cur_sur_top_bou_tinsId_np = np.array(cur_sur_top_bou_tinsId_np)
        cur_sur_OutcropBou_pds = pd.DataFrame({
            'xx1':list(cur_sur_top_bou_tinsNode_np[:,0].flatten()),
            'yy1':list(cur_sur_top_bou_tinsNode_np[:,1].flatten()),
            'zz1':list(cur_sur_top_bou_tinsNode_np[:,2].flatten()),
            'xx2':list(cur_sur_top_bou_tinsNode_np[:,3].flatten()),
            'yy2':list(cur_sur_top_bou_tinsNode_np[:,4].flatten()),
            'zz2':list(cur_sur_top_bou_tinsNode_np[:,5].flatten()),
            'xx3':list(cur_sur_top_bou_tinsNode_np[:,6].flatten()),
            'yy3':list(cur_sur_top_bou_tinsNode_np[:,7].flatten()),
            'zz3':list(cur_sur_top_bou_tinsNode_np[:,8].flatten())
        })
        same11 = cur_sur_OutcropBou_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'last')#保留最后一个，，其他的标记True
        same22 = cur_sur_OutcropBou_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'first')#除第一个其他标记True
        same_all2 = np.array(list(zip(list(same11),list(same22))))
        same_all_result2 = same_all2.any(axis = 1)#获取所有的True，后面提取这些重复的点
        same_all_result2 = np.hstack((same_all_result2[np.arange(cur_sur_tinsNode_count)],np.zeros(same_all_result2.shape[0]-cur_sur_tinsNode_count,dtype=bool)))
        trgl_result2 = cur_sur_top_bou_tinsId_np[same_all_result2]
        OutCropBou = MultiPolygon()
        Result_top_bou_Surface[lay_name]['xyz9'] = list(cur_sur_top_bou_tinsNode_np[same_all_result2])
        Result_top_bou_Surface[lay_name]['id3'] = list(cur_sur_top_bou_tinsId_np[same_all_result2])
        if trgl_result2.size>0:
            trgl_result2_p1 = trgl_result2[:,0]
            trgl_result2_p2 = trgl_result2[:,1]
            trgl_result2_p3 = trgl_result2[:,2]
            nodes_col_1 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p1-1][:,:3]
            nodes_col_2 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p2-1][:,:3]
            nodes_col_3 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p3-1][:,:3]
            Line1 = np.hstack((nodes_col_1,nodes_col_2))
            Line2 = np.hstack((nodes_col_1,nodes_col_3))
            Line3 = np.hstack((nodes_col_2,nodes_col_3))
            LineAll = np.vstack((Line1,Line2))
            LineAll = np.vstack((LineAll,Line3))
            Line_init = deepcopy(LineAll)
            LineAll.sort(axis=1)
            OutCrop_Bou_pd = pd.DataFrame({
                'xx1':list(LineAll[:,0].flatten()),
                'yy1':list(LineAll[:,1].flatten()),
                'zz1':list(LineAll[:,2].flatten()),
                'xx2':list(LineAll[:,3].flatten()),
                'yy2':list(LineAll[:,4].flatten()),
                'zz2':list(LineAll[:,5].flatten()),
            })
            same11 = OutCrop_Bou_pd.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2'],keep = 'last')#保留最后一个，，其他的标记True
            same22 = OutCrop_Bou_pd.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2'],keep = 'first')#除第一个其他标记True
            same_all2 = np.array(list(zip(list(same11),list(same22))))
            same_all_result2 = same_all2.any(axis = 1)#获取所有的True，后面提取这些重复的点
            trgl_same_reverse2 = np.logical_not(same_all_result2)
            Line_part = Line_init[trgl_same_reverse2]
            if Line_part.size>0:
                if Line_part.ndim ==2:
                    for l_i in Line_part:
                        xx1,yy1,zz1,xx2,yy2,zz2 = l_i
                        OutCropBou = OutCropBou.union(LineString([[xx1,yy1,zz1],[xx2,yy2,zz2]]))
            if OutCropBou:
                f_out2 = open(FolderPath+'/TopBotSurface/'+lay_name+'_OutCropBou.dat','w')
                
                if OutCropBou.geom_type == "LineString" or OutCropBou.geom_type == "LinearRing":
                    f_out2.write('line 1 coords\n')
                    coords_bou = np.array(OutCropBou.wkt[OutCropBou.wkt.find('(')+1:-1].replace(', ',',').replace(' ',',').split(','))
                    coords_bou = coords_bou.reshape((-1,3))
                    for bou_coord in coords_bou:
                        f_out2.write(bou_coord[0]+','+bou_coord[1]+','+bou_coord[2]+'\n')
                elif OutCropBou.geom_type == "MultiLineString":
                    for boui in range(len(OutCropBou.geoms)):
                        bou_temp = OutCropBou.geoms[boui]
                        f_out2.write('line '+str(boui)+' coords\n')
                        coords_bou = np.array(bou_temp.wkt[bou_temp.wkt.find('(')+1:-1].replace(', ',',').replace(' ',',').split(','))
                        coords_bou = coords_bou.reshape((-1,3))
                        for bou_coord in coords_bou:
                            f_out2.write(bou_coord[0]+','+bou_coord[1]+','+bou_coord[2]+'\n')
                f_out2.close()

        print(lay_name+' Outcrop Boundary extract OK...')
        #Get Eroded Bou of the layer
        
        Result_Erod_bou_Surface ={}
        Result_Erod_bou_Surface[lay_name]={}
        Result_Erod_bou_Surface[lay_name]['nodes']=[]
        Result_Erod_bou_Surface[lay_name]['tins']=[]
        Result_Erod_bou_Surface[lay_name]['xyz9']=[]
        Result_Erod_bou_Surface[lay_name]['id3']=[]
        cur_sur_Erod_bou_tinsNode_np = deepcopy(Surface_All[lay_name]['xyz9'])
        cur_sur_Erod_bou_tinsId_np = deepcopy(Surface_All[lay_name]['id3'])
        
        for erii in range(after_Eroded_count):
            cur_sur_Erod_bou_tinsNode_np.extend(Result_Eroded_Surface[erii]['xyz9'])
            cur_sur_Erod_bou_tinsId_np.extend(Result_Eroded_Surface[erii]['id3'])

        cur_sur_Erod_bou_tinsNode_np = np.array(cur_sur_Erod_bou_tinsNode_np)
        cur_sur_Erod_bou_tinsId_np = np.array(cur_sur_Erod_bou_tinsId_np)
        cur_sur_ErodBou_pds = pd.DataFrame({
            'xx1':list(cur_sur_Erod_bou_tinsNode_np[:,0].flatten()),
            'yy1':list(cur_sur_Erod_bou_tinsNode_np[:,1].flatten()),
            'zz1':list(cur_sur_Erod_bou_tinsNode_np[:,2].flatten()),
            'xx2':list(cur_sur_Erod_bou_tinsNode_np[:,3].flatten()),
            'yy2':list(cur_sur_Erod_bou_tinsNode_np[:,4].flatten()),
            'zz2':list(cur_sur_Erod_bou_tinsNode_np[:,5].flatten()),
            'xx3':list(cur_sur_Erod_bou_tinsNode_np[:,6].flatten()),
            'yy3':list(cur_sur_Erod_bou_tinsNode_np[:,7].flatten()),
            'zz3':list(cur_sur_Erod_bou_tinsNode_np[:,8].flatten())
        })
        same11 = cur_sur_ErodBou_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'last')#保留最后一个，，其他的标记True
        same22 = cur_sur_ErodBou_pds.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2','xx3','yy3','zz3'],keep = 'first')#除第一个其他标记True
        same_all2 = np.array(list(zip(list(same11),list(same22))))
        same_all_result2 = same_all2.any(axis = 1)#获取所有的True，后面提取这些重复的点
        same_all_result2 = np.hstack((same_all_result2[np.arange(cur_sur_tinsNode_count)],np.zeros(same_all_result2.shape[0]-cur_sur_tinsNode_count,dtype=bool)))
        trgl_result2 = cur_sur_Erod_bou_tinsId_np[same_all_result2]
        if trgl_result2.size>0:
            ErodBou = MultiPolygon()
            trgl_result2_p1 = trgl_result2[:,0]
            trgl_result2_p2 = trgl_result2[:,1]
            trgl_result2_p3 = trgl_result2[:,2]
            nodes_col_1 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p1-1][:,:3]
            nodes_col_2 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p2-1][:,:3]
            nodes_col_3 = np.array(Surface_All[lay_name]['nodes'])[trgl_result2_p3-1][:,:3]
            Line1 = np.hstack((nodes_col_1,nodes_col_2))
            Line2 = np.hstack((nodes_col_1,nodes_col_3))
            Line3 = np.hstack((nodes_col_2,nodes_col_3))
            LineAll = np.vstack((Line1,Line2))
            LineAll = np.vstack((LineAll,Line3))
            Line_init = deepcopy(LineAll)
            LineAll.sort(axis=1)
            Erod_Bou_pd = pd.DataFrame({
                'xx1':list(LineAll[:,0].flatten()),
                'yy1':list(LineAll[:,1].flatten()),
                'zz1':list(LineAll[:,2].flatten()),
                'xx2':list(LineAll[:,3].flatten()),
                'yy2':list(LineAll[:,4].flatten()),
                'zz2':list(LineAll[:,5].flatten()),
            })
            same11 = Erod_Bou_pd.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2'],keep = 'last')#保留最后一个，，其他的标记True
            same22 = Erod_Bou_pd.duplicated(subset = ['xx1','yy1','zz1','xx2','yy2','zz2'],keep = 'first')#除第一个其他标记True
            same_all2 = np.array(list(zip(list(same11),list(same22))))
            if same_all2.size>0:
                same_all_result2 = same_all2.any(axis = 1)#获取所有的True，后面提取这些重复的点
                trgl_same_reverse2 = np.logical_not(same_all_result2)
                Line_part = Line_init[trgl_same_reverse2]
                if Line_part.size>0:
                    if Line_part.ndim ==2:
                        for l_i in Line_part:
                            xx1,yy1,zz1,xx2,yy2,zz2 = l_i
                            ErodBou = ErodBou.union(LineString([[xx1,yy1,zz1],[xx2,yy2,zz2]]))
                if ErodBou:
                    f_out2 = open(FolderPath+'/TopBotSurface/'+lay_name+'_ErodBou.dat','w')
                    
                    if ErodBou.geom_type == "LineString" or ErodBou.geom_type == "LinearRing":
                        f_out2.write('line 1 coords\n')
                        coords_bou = np.array(ErodBou.wkt[ErodBou.wkt.find('(')+1:-1].replace(', ',',').replace(' ',',').split(','))
                        coords_bou = coords_bou.reshape((-1,3))
                        for bou_coord in coords_bou:
                            f_out2.write(bou_coord[0]+','+bou_coord[1]+','+bou_coord[2]+'\n')
                    elif ErodBou.geom_type == "MultiLineString":
                        for boui in range(len(ErodBou.geoms)):
                            bou_temp = ErodBou.geoms[boui]
                            f_out2.write('line '+str(boui)+' coords\n')
                            coords_bou = np.array(bou_temp.wkt[bou_temp.wkt.find('(')+1:-1].replace(', ',',').replace(' ',',').split(','))
                            coords_bou = coords_bou.reshape((-1,3))
                            for bou_coord in coords_bou:
                                f_out2.write(bou_coord[0]+','+bou_coord[1]+','+bou_coord[2]+'\n')
                    f_out2.close()

        print(lay_name+' Eroded Boundary extract OK...')
    print('......Run completed.....')
    

        




    