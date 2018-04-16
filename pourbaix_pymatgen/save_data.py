def save_data(data,file_name):
	"""
	Saves data into a json format
	Args:
		data: data to be saved
		file_name: Desired filename including extension (.json)
	"""
	#| -  - save_data
	import json
	# Storing and Retrieving Matrix data
	json_dumps = json.dumps(data)
	f = open(file_name,'w')
	f.write(json_dumps); f.close()
	print 'Saving file: '+str(file_name)
	#__|
