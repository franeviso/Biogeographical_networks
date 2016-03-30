# execfile("Y:\Gython Script.py")
import re
from java.awt import Color
Spectral = {
    3: ['#99D594', '#FFFFBF', '#FC8D59'],
    4: ["#2B83BA", "#ABDDA4", "#FDAE61", "#D7191C"],
    5: ["#2B83BA", "#ABDDA4", "#FFFFBF", "#FDAE61", "#D7191C"],
    6: ["#3288BD", "#99D594", "#E6F598", "#FEE08B", "#FC8D59", "#D53E4F"],
    7: ['#3288BD', '#99D594', '#E6F598', '#FFFFBF', '#FEE08B', '#FC8D59', '#D53E4F'],
    8: ['#3288BD', '#66C2A5', '#ABDDA4', '#E6F598', '#FEE08B', '#FDAE61', '#F46D43', '#D53E4F'],
    9: ['#3288BD', '#66C2A5', '#ABDDA4', '#E6F598', '#FFFFBF', '#FEE08B', '#FDAE61', '#F46D43', '#D53E4F'],
    10: ["#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"],
    11: ['#5E4FA2', '#3288BD', '#66C2A5', '#ABDDA4', '#E6F598', '#FFFFBF', '#FEE08B', '#FDAE61', '#F46D43', '#D53E4F', '#9E0142'],
	33: ['#5E4FA2', '#5060AA', '#4272B2', '#3484BB', '#3E96B7', '#4FA8AF', '#5FBAA8', '#72C7A4', '#88CFA4', '#9ED7A4', '#B2E0A2', '#C4E79E', '#D7EF9B', '#E7F59A', '#EFF8A6', '#F7FBB2', '#FFFFBF', '#FEF5AE', '#FEEB9E', '#FEE18E', '#FDD380', '#FDC373', '#FDB466', '#FBA15B', '#F88D52', '#F57948', '#F06744', '#E65848', '#DC494C', '#D13A4E', '#C0274A', '#AF1446', '#9E0142'],
	38: ['#5E4FA2', '#525EA9', '#466DB0', '#3A7DB7', '#368CBB', '#449CB4', '#52ACAE', '#60BBA7', '#71C6A4', '#83CDA4', '#96D4A4', '#A9DCA4', '#B9E2A1', '#C9E99D', '#D9EF9A', '#E7F59A', '#EEF8A4', '#F4FAAF', '#FBFDB9', '#FEFAB7', '#FEF2A9', '#FEEA9B', '#FEE18D', '#FDD581', '#FDC776', '#FDBA6B', '#FCAC60', '#FA9A58', '#F7894F', '#F57747', '#F06744', '#E85B47', '#DF4E4A', '#D7414E', '#CA324C', '#BB2149', '#AC1145', '#9E0142'],
    41: ['#5E4FA2', '#535DA8', '#486BAF', '#3C79B6', '#3288BD', '#3E96B7', '#4CA5B1', '#59B3AB', '#66C2A5', '#77C8A4', '#88CFA4', '#99D6A4', '#ABDDA4', '#B9E3A1', '#C8E99E', '#D7EF9B', '#E6F598', '#ECF7A1', '#F2FAAB', '#F8FCB5', '#FFFFBF', '#FEF7B2', '#FEEFA4', '#FEE798', '#FEE08B', '#FDD380', '#FDC776', '#FDBA6B', '#FDAE61', '#FA9D59', '#F88D52', '#F67D4A', '#F46D43', '#EC6146', '#E45549', '#DC494C', '#D53E4F', '#C72E4B', '#B91F48', '#AB1045', '#9E0142']
}


cluster_type = 'map_equation_rw_moduleID'
cluster_type = 'modularity_conv_moduleID'
cluster_type = 'map_equation_conv_moduleID'


# runLayout(ForceAtlas2, iters = 3000)

#num_clusters = 10
#cluster_numbers = [1, 3, 4, 5]
# cluster_numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
cluster_numbers = []
for node in list(g.nodes):
	cluster_num = int(getattr(node, cluster_type))
	if (cluster_num not in cluster_numbers):
	
		cluster_numbers.append(cluster_num)
cluster_numbers.sort()
print cluster_numbers
num_clusters = len(cluster_numbers)

colour_list = Spectral[num_clusters]

for i in range(num_clusters):
#for i in range(num_clusters):
	mySubGraph = g.filter(eval(cluster_type) == cluster_numbers[i])
	print i
	print cluster_numbers[i]
	print colour_list[i]
	for node in list(mySubGraph.nodes):
		node.color = Color.decode(colour_list[i])
		
# print rgb
# node.color = Color(rgb[0], rgb[1], rgb[2])
# rgb = [int(s) for s in re.split('[,]|[(]|[)]', colour_list[i]) if s.isdigit()]