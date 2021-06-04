import pickle, requests, json, mne, ctypes
import matplotlib
from mne.viz import circular_layout, plot_connectivity_circle
import numpy as np
import matplotlib.pyplot as plt
import PySimpleGUI as sg
import mne.viz.utils as utils
from functools import partial

def export(item, subPath, name):
    mainPath = 'C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE'
    file = open(mainPath + subPath + '/%s.pkl' % name, 'wb')
    print('saving to' + str(file))
    pickle.dump(item, file)
    file.close()

def get(subPath, name):
    mainPath = 'C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE'
    file = open(mainPath + subPath + '/%s.pkl' % name, 'rb')
    item = pickle.load(file)
    file.close()
    return item

class formatData():
    def __init__(self):
        self.sym2ensg = "/xrefs/symbol/canis_lupus_familiaris/"
        self.orthology = "/homology/symbol/canis_lupus_familiaris/"
        self.orthologyId = "/homology/id/"
        self.translate = "/map/translation/%s/1..100000?"
        self.server = "http://rest.ensembl.org"
    
    def buildCorrelationMap(listA, listB, length):
        print('finding correlations for genes in slice with other slices %s' % listA)
        cor = []
        string = 'https://string-db.org'
        interactions = '/api/tsv/network?identifiers='
        for b, y in enumerate(listB):
            if listA == y[0]:
                continue
            total = 0.0
            correlations = requests.get(string+interactions+listA+y[0]+'&species=9612')
            correlations = correlations.text
            correlations = correlations.split('\n')
            for c, line in enumerate(correlations):
                line = line.split('\t')
                if len(line) > 5:
                    if line[5] != 'score':
                        # print(line[2],line[3],line[5])
                        if line[2] in listA or line[3] in listA:
                            total += float(line[5])
            cor.append(float(total)/c)
        if len(cor) < length:
            for m in range(length-len(cor)):
                lst = []
                for i in range(length):
                    lst.append(0.0)
                cor.append(lst)
        print(cor)
        return cor

    def mainCorrelationMap(geneList, slices, location):
        length = len(slices)
        corMap = []
        labelNames = []
        labelColors = []
        for i in range(length-1):
            label = slices[i][2]
            #labelNames.append(str(label[0]) + '@' + location[i])
            labelNames.append(str(label[0])+str(i))
            a = label[1]
            labelColors.append((0.75,0.3,0.1,a))
        for i in range(length-1):
            corMap.append(formatData.buildCorrelationMap(geneList[i][0], geneList, length-1))
        #for i in range(length-1):    
        #    if i <5:
        #        corMap.append(formatData.buildCorrelationMap(geneList[i][0], geneList, length-1))
        #    else:
        #        lst = []
        #        for m in range(length-1):
        #            lst.append(0.0)
        #        corMap.append(lst)
        print(corMap)
        return corMap, labelNames, labelColors

    def prepareHomologyGraph(self, labelNames, node, geneList, sliceCorrelations):
        genesInSlice = geneList[node]
        descriptions = genesInSlice[1].split('\n')
        genesInSlice = genesInSlice[0].split('%0d')
        homologs = {}
        filterOrthologs = input('type True to filter orthologs or False to filter paralogs')
        if json.loads(filterOrthologs.lower()) == True:
            typeFilter = 'ENSEMBL_ORTHOLOGUES'
        else:
            typeFilter = 'ENSEMBL_PARALOGUES'
        print(typeFilter)
        for i, a in enumerate(genesInSlice):
            b = descriptions[i]
            ENSG = requests.get(self.server+self.sym2ensg+a, headers={ "Content-Type" : "application/json"})
            if ENSG.ok:
                ENSG = ENSG.json()
                #print(ENSG)
                print(a, b)
                for i in ENSG:
                    print(i["id"])
                    homology = requests.get(self.server+self.orthologyId+i["id"], headers={ "Content-Type" : "application/json"})
                    if homology.ok:
                        homology = homology.json()
                        #print(str(homology))
                        homologs[a] = []
                        homologs[a].append(['id', 'protein_id', 'perc_id', 'per_pos', 'cigar_line', 'align_seq'])
                        for x in homology:
                            if x != 'error':
                                homologies = homology['data'][0]['homologies']
                                for A, data in enumerate(homologies):
                                    #print(data)
                                    #print('\n\nmethod_link_type:\t', data['method_link_type'])
                                    #print(typeFilter)
                                    if data['method_link_type'] != typeFilter:
                                        print('\n\nmethod_link_type:\t', data['method_link_type'])
                                        for mapping in homologies[A]:
                                            #print(mapping)
                                            if mapping == 'source':
                                                if float(homologies[A][mapping]['perc_id']) > 40 and float(homologies[A][mapping]['perc_pos']) > 50:
                                                    print(mapping)
                                                    print('id:\t', homologies[A][mapping]['id'])
                                                    print('protein_id:\t', homologies[A][mapping]['protein_id'])
                                                    print("perc_id:\t", homologies[A][mapping]['perc_id'])
                                                    print("perc_pos:\t", homologies[A][mapping]['perc_pos'])
                                                    print('cigar_line:\t', homologies[A][mapping]['cigar_line'])
                                                    print('align_seq:\t', homologies[A][mapping]['align_seq'])
                                                    homologs[a].append([mapping, homologies[A][mapping]['id'], homologies[A][mapping]['protein_id'], homologies[A][mapping]['perc_id'], homologies[A][mapping]['perc_pos'], homologies[A][mapping]['cigar_line'], homologies[A][mapping]['align_seq']])
                                            elif mapping == 'target':
                                                if float(homologies[A][mapping]['perc_id']) > 40 and float(homologies[A][mapping]['perc_pos']) > 50:
                                                    homologs[a].append([mapping, homologies[A][mapping]['id'], homologies[A][mapping]['protein_id'], homologies[A][mapping]['perc_id'], homologies[A][mapping]['perc_pos'], homologies[A][mapping]['cigar_line'], homologies[A][mapping]['align_seq']])
                                                    
                                break
        print(genesInSlice)
        index = input('type Gene symbol to print homologies')
        print(homologs[index])
        for i, data in enumerate(homologs[index]):
            if i > 1:
                print(data, self.server+self.translate % data[2])
                transcriptCoords = requests.get(self.server+self.translate % data[2], headers={ "Content-Type" : "application/json"}) 
                if transcriptCoords.ok:
                    transcriptCoords = transcriptCoords.json()
                    print(str(transcriptCoords))
        thisLabel = labelNames[node]
        print(sliceCorrelations)

    def prepareCorrelationGraph(self, labelNames, node, geneList, sliceCorrelations):
        print('pickling')
        genesInSlice = geneList[node]
        print(genesInSlice)
        string = 'https://string-db.org'
        interactions = '/api/tsv/network?identifiers='
        for i, percent in enumerate(sliceCorrelations):
            print(percent)
            if percent > 0:
                correlations = requests.get(string+interactions+geneList[i][0]+genesInSlice[0]+'&species=9612')
                correlations = correlations.text
                print(correlations)
            #correlations = correlations.split('\n')
        descriptions = genesInSlice[1].split('\n')
        genesInSlice = genesInSlice[0].split('%0d')
        for i, gene in enumerate(genesInSlice):
            print(gene)
            print(descriptions[i])
        #subPath = '/DogNames/Masha5'
        #correlations = {'genesInSlice':  open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE'+subPath+'/genesInSlice.pkl', 'wb'), 'geneList':  open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE'+subPath+'/geneList.pkl', 'wb'), 'sliceCorrelations': open('C:/Users/webaccus/OneDrive - Agilent Technologies/Desktop/DGAE'+subPath+'/sliceCorrelations.pkl', 'wb')}
        #pickle.dump(genesInSlice, correlations['genesInSlice'])
        #correlations['genesInSlice'].close()
        #pickle.dump(geneList, correlations['geneList'])
        #correlations['geneList'].close()
        #pickle.dump(sliceCorrelations, correlations['sliceCorrelations'])
        #correlations['sliceCorrelations'].close()

class displayGraphics():
    def __init__():
        import matplotlib.pyplot as plt

    def grabData(event, fig=None, axes=None, indices=None,
                                         n_nodes=0, node_angles=None,
                                         ylim=[9, 10], labelNames=None, geneList=None, corMap=None):
        if event.inaxes != axes:
            return

        #event buttons: 1 = left click, 2 = double click, 3 = right click
        if event.button == 3:  # right click
            # click must be near node radius
            if not ylim[0] <= event.ydata <= ylim[1]:
                return
            node_angles = node_angles % (np.pi * 2)
            node = np.argmin(np.abs(event.xdata - node_angles))
            return formatData().prepareHomologyGraph(labelNames, node, geneList, corMap[node])
        elif event.button == 1:  # left click
            # click must be near node radius
            if not ylim[0] <= event.ydata <= ylim[1]:
                return
            node_angles = node_angles % (np.pi * 2)
            node = np.argmin(np.abs(event.xdata - node_angles))
            return formatData().prepareCorrelationGraph(labelNames, node, geneList, corMap[node])

    def buildConnectivityMap(slices, length, geneList, corMap, labelNames, labelColors):
        node_order = list()
        node_order.extend(labelNames)
        node_angles = circular_layout(labelNames, node_order, start_pos=90, 
                                  group_boundaries=None)
        label_names = np.array(labelNames)
        label_colors = np.array(labelColors)
        cor_map = np.array(corMap)
        #instance = x+1
        print('Almost Done! Building connectivity circle and UI')
        fig, axes, indices, n_nodes, node_angles, ID, plt = plot_connectivity_circle(cor_map, label_names, n_lines=300, node_angles=node_angles, 
                             node_colors=label_colors,
                             title='Chromosome', show=False)
        print(fig, axes)
        callback = partial(displayGraphics.grabData, fig=fig,
                               axes=axes, indices=indices, n_nodes=n_nodes,
                               node_angles=node_angles, labelNames=labelNames, geneList=geneList, corMap=corMap)
        fig.canvas.mpl_connect('button_press_event', callback)
        utils.plt_show(True)
        #plt.show()

    def plotSubChromosome():
        print('nope')

#dogName = 'Masha'
#chromosomeNumber = 1
#filename = '%s/%d/' % (dogName, chromosomeNumber)
# Layout the design of the GUI
def showLayout():
    chromosome = b'iVBORw0KGgoAAAANSUhEUgAAAB4AAAA3CAYAAAAbk/pcAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAA0ISURBVFhHfZhriF3Xdcf3Pq/7mLnzliVZ0ki2ZcV20jSkqV1jNzh287JDwIE8oIQUQmmhycdS8iHkYz/0Q+mHltLS0NKGEhJapcWFvFwH4iTUMU6TYBQpkvWwNJLmpdGdua/z2P3997nnzpnR0MU9Z++z93rsvfZaa6917bPPPh0b43Jr87CzeKoRJ+2HjDGPWxuERZ5dTEeDV7obN0wQ5MYBRRFa5wIXRSP1g4rWOVs4F4Hngs78/S5K2r8VBOFvOuOE85bL83PrN9+6UdHap59+Mmg2o3DhyJnnnDWfROi7rbGHQbbG2CGyfobEv05H2693N673EGAsoovCBWnqctHSz8EPp+ePxkmz86i19lMs/OOMTcMDli7jdZX2m1k6+trmrUvd8PTpU9HC0TNPgPyHoHwEtg+CMAfBDO0C48esNaeDML7sTHE1T4cOJkGWFXkch+zUFOCFUTJl2jNLjyHw8zyfYewkPGZp4WPgZ45AdgRFNrt3Vl8Nmq1p7ewFJt5He0hM6APi7wFi8xTM3t+aWjxO328hDANASI79s7W5+2bA+W0W+jHGEOLx6tBg/B2MPr949PR80J45JMbP8ty/K6yi8TLUafJw7uGjOmdBECCG0xMgzKCREzSP8rk8HuWZLB7wPGfQ6CNhFD0WQHEGwhMMtsvJiqBO6MePIgGmEovCQqwNw5FRaR4e2IVBIwfRCibjTVT1OzwhZ4E5TpArgjqhH0eQKRAAjHKkR72eTRkLsywXQoTGpRlgP63aSZ/FutN6cR6Wc62Q67BLyPsWnVu4DuebhP2+S+PYJlh4gZHBp7jOvt8SVQl1oWon/ZBDOsxJ2Sk+NLIP9hLy/hXMbxVFYPPcFM2mwaZthsbgkbs8G10G721R3itUrcD3tdlpXrI2i5loskIQ1AkN/mxw/vxqUeiQdaQWtRdYdEZscCYdoANntONr+2jHbcXfEyeyy5q0CrE2VMJtnuuj/s5mGUDAUBRxBbsvIApNr7suY1sF76Kn2AOlwDF/fXiH2ClH6kLV323ZyTVeq1vrN4elGwVMFEEuqU6nFZhQ3u/cGriofA/tpF8CXuDcUILXmGC1JZL0WCdUyw5X6N6ltZyrw58JkxmoAYuItHlwpPZiB5xf76WtBE7Gcpp1BLt1BuQeHqlE3CucqZT5LIq4AZCFOxX9vlTcMIwVEXcDOPzCEa9hnXa3rfoGLZmbgSvyHzKQ7QosEaudlwvScXCg6FTGlKaFaTSmbJJkXBZpOBq5TKaFZ0Fj/CbqPHb5esD33ZtBr7u2xSRu4EYlEj0vrCRQy+cJdNGSerNMWol5+kWe26AorHQe4tOQ21neR+q0avcBcowE35Eaf8bHhka1SI/vCSoG5iSvJa4/P8edhFkSM3Ocy3ApowdUjhy7AMnpOu0+UKRbL/LRRVkmffc9XnfGQgD/qsNx9nl49tBywhn7AXZqRCuvkqW3Zg5JDUeZ0kVxAHg6ZJj/HWxvruEZ+EKRv8YA1mjH6t6/YqKbNe+Jk9ZDch2dxM5OVogWN8LKrYsb04u03OXm1D7acesPfAXtvry+up4GrZYpttYuXGfwDSYJFPuFCnz/o8h5odFeSAaDglBpPS3Cgqm5YzGqfwKcZ6BtHUAL2AGvXyP+u6LFOOIoTcMc8zzLxE/vFSrwq9Y5P9mann+XjGl2tgBK2ihuLbMAUh2L4ANpAXee13dWLp9bEW14/Pgyq3Zxr7t6s91ZmgXvJKjcrRVBtRD//QACDndml87t3N3eGI2K0fzhE0eiOPkz5j7Is1Ai3kPb5XU2z9KvjvprXedacRBFLsI2yRRDjjr/Nnj/BhIuVgmrt97yN4hQfe7gdP6+EwtkpV9h7P08tXRnD+0ADZ1Fxd8b9e9u4Iq6v9Pw1KnjMhA2AkoQ78RxY40PEjT7Hgg1KgZEI/M6zzdd4c7iz5emZhZPstMvMf0840p3Ep4xTIQqEfwOzz/mWfZab+v6Dt7ARBGQZRIbwNCtM+z3iqQ51Q3CaJWFKDt8jIdYbl7GKP+dqPjf7PYa+dW7sKlPg/NZ5pZ4yM0FYjUB+Sz47p+LPP3J3bXLq5KhK1RxHcEnUbUCrb8h7bC/lSWtKVQSXoW5BP8PxF8rivxVMbNB9BQtQu3v02pxENYFakdGl84vofsrioKXEbouF1TIYLEgp85+4AO/G+V5Rjqj0EP2Fij426A9e79LGi2qCns7S3sFefMpnODDMP0ctCxI/AWVWgWWxN1t067iJV/O0sHZ7sZlLgVLshDg7wYZKZdKEthnnnkq5Ix1myqEoXLWk5ZJehj2sjCcThaPnnmBpf4J80+CUvPTSuhEOImAfSnLRn++sfLzC6LFcFO0o1sQJPyQAJvncTGunQqsuvD1D8kiMpwNk0bRbC/MNNuzX4bgRZjKapVFSsq42SP0v2i+nufpy6PB9s0szZLe9t0BBQ4Xti4sSo5C+wsQnjr8+LgpSxGfH8u6bauzFLWml+5Lmu0/BfPTjCsGHygU69jk659kB1yxrxNEcbHWi43W9BfandlLcdzcGg22OAJ/9kgoOAUSiocfPqn6RwvyUtudxbg5vfAIZ/4H8EUoEUtJ1T3ghSoafRXib6O0c2hUsfzj8BHd+2hn8JCrYdRYHw12dJQSoizTyarp5FbZY3t6LmzNLD3I3O+xhj8Cbxzw6+C/5dfnkfwv7PIb9LcwvN/grnwBuk+A8w7GpKFHEITTxudH/e2NsiaQbHby4IPLND5jc52FY0sE++fo/zEIZ0ohXkNj8N9Kbc6z6L8ll/77LO13qU0/DLsvMvk8J4VfV2fv/RvPCO4kzemLw97mtnarXXNJFCRued5oz4SsTOr5EMjv3ivUMwF8loLQ4l9vXzv/N3duXxrEjSkKAvMxMJ5EKLushE5oF/n8LLxfJF4gNCBkmiw8duxkgd9Es0snDrHbz4P5KZ5xJBJUTNRaqgn3d6tv//wvkySOh0ObDnsbA+ziBgtW0fYwTw1KWmkB9Q6wnXP97spKHEdxQK5U/hURRB8FixrZn80Y6kJdn+eH5Hzfj6JmmKYuq2gHO3dewzy/D9IFnjHsodXAR8g6Pqe6S7T0Tb5w5CGqfvscSJyrQqDgHsJX2e2PR4PuFQyR7MMRCkxOnh327q6MOLL/BO0lId5Lq9Zv6L1zhx9+b0lLJCMpfycI1MiuM0YC9hCOaH7Krt64u7Ei4xIQXhWNFGJJHYe9t3EvEglVlfuF+lYdiv/gWU8rI2PscQZ1yyi0CYmnDm6d0V+lw/6l8lpTpon7Q6vaWBfMzp1bGJ67wRS5Wx08+riVxZsnRBskjY7Wwt3rswfgnpWyUXOF1+3N1bd72irWSZQzxNxUuTVTgf5mEjpJgk8a67Rj8GMdyN9JjjYTNKcWpxlBsD0kJJZfIU1aFqZkcJO+jSLyJMySW0Y7pYqKdauB42unLoz/v9pJfwDMRnHz8SCIE6KMbUptQioR9wpnLVzfLq/XTtvbGZNU53FYxHFVOwWqKhTz7+Gx2yf8smslmscYUOQaT5SI1c7VMtdRcKAl4GSO9NYlyZRtNgeqnfwfbb52soYr03TqtLsC9QhsxOuwsnMOvP4fyEQY/bLlW5lGooSBqEOFiAPHVGpZxDmThaFzdi0KuYwS+wltyVdt9WiT5rh2qot9T/ri8ScEYmAeoDurbzHDd4UQZBSLfKM0XzuJcJ5paXBCuxf8t27DBQWL8ZVXIVVEe+A+Ro7NLS23ySD8wuRWg0EqX8airWm058TnELTcaPuhvgC/rVbpx+XHGNSthFcEnAt1URg3fAkq1Xa7yp2IP+RoGktas21a3UwH1E67PMfHKFewvuguJwX7hQp8/xj2/ICyh14vTyWy05GqpWhFsWgWnkQm2z6AFih5apG0wwCzJG/2ifd4cr9QgVaqujdYHgwybnMXzs35yzRW7TQcUrP67KOqjetQ8av4O+Xbb8r/9FcEH3UEgb6r1q9Uf7Z2ksQF7TYEaStiEfpbIVQehdAOlLJ+YC/t7rfv6x/28wFXGsWyvcKgysgx1AkmrQqvniwaG9a/AbnUjZp1nyst3gGLYzuQlkfACo3T31s/CLa3bmnV1Df+T7Qx7Cfw7UUkXtMZKT/L5NCk38T78bG5ax7nYNoxOIpB86N0NLhIwUaMKwpViKhcVaKI6uC/9Qcbl31xAaEmy3Am5d5KimUeLIBc+goW+wofl0RQQp2XldDXudv+Y+PW5QHp7alo5+7aanN6aYulzzOpSKZSRJWOLF5X3bdQ0UvdzZtvjobyXcO5DrksEkUu1hRQ8G0XzdYM2WZ4Ew08Aq2MUAKVkarw+wWf3+De/tawf3sYLi+f1Hm5wc7WldbU/C+xIv8PLIvYhOMbdP8iG/W/vnnr4oV0lPvLP+AuiOMWMZrVQYsrcdokA4ONHsn8JdKbVxC5ozNAC+L3XZ5/WL9xHqHrwzhuBf8HcMzEaJ5SIxIAAAAASUVORK5CYII='
    layout = [[sg.Text('Please click a new button', auto_size_text=True)]]
    x = 1
    layoutGrid = []
    end = 11
    for i in range(1, 40):
        if i == end:
            end += 10
            layout.append(layoutGrid)
            layoutGrid = []
        layoutGrid.append(sg.Button('%d' % i, image_data=chromosome,
              button_color=('black',sg.theme_background_color()), border_width=0, key='%d' % i))
        if i == 39:
            layout.append(layoutGrid)
    return layout

def viewLoop(window, dogName):
    from mne.viz import circular_layout, plot_connectivity_circle
    import mne.viz.utils as utils
    import numpy as np

    subPath = '/DogNames'

    # Lookup dictionary that maps button to function to call
    dispatch_dictionary = {'1':str(1), '2':str(2), '3':str(3), '4':str(4), '5':str(5), '6':str(6), '7':str(7), '8':str(8), '9':str(9), '10':str(10), '11':str(11), '12':str(12), '13':str(13), '14':str(14), '15':str(15), '16':str(16), '17':str(17), '18':str(18), '19':str(19)}
    
    # Event loop. Read buttons, make callbacks
    while True:
        # Read the Window
        event, value = window.read()
        if event in ('Quit', sg.WIN_CLOSED):
            break
        # Lookup event in function dictionary
        if event in dispatch_dictionary:
            chromosome = dispatch_dictionary[event]   # get function from dispatch dictionary
            file = '/'+dogName+'/'+str(chromosome)+'/'
            corMap, labelNames, labelColors = formatData.mainCorrelationMap(get(subPath, file + 'geneList'), get(subPath, file + 'alignedSlices'), get(subPath, file + 'location'))
            print(len(labelNames), len(corMap))
            return displayGraphics.buildConnectivityMap(get(subPath, file + 'alignedSlices'), get(subPath, file + 'length'), get(subPath, file + 'geneList'), corMap, labelNames, labelColors)
        else:
            print('Event {} not in dispatch dictionary'.format(event))
