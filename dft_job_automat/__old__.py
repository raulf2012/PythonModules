#| - job_analysis

#| - DEPRECATED METHODS
# def max_force(self, path):
#     """
#     """
#     #| - max_force
#     with open(path + "/simulation/qn.log", "r") as file:
#         max_f = file.readlines()[-1].split(" ")[-1]
#         return(float(max_f))
#     #__|

# def job_id_names(self):
#     """
#     """
#     #| - job_id_names
#     try:
#         f = open("job_id_names", "r")
#         job_id_data = pickle.load(f)
#         f.close()
#
#         id_list = job_id_data[0]
#         name_list = job_id_data[1]
#
#         self.data_frame["job_id"] = id_list
#         self.data_frame["job_name"] = name_list
#     except:
#         print("Couldn't parse job_id_names file")
#
#     #__|
#__|

#__|
