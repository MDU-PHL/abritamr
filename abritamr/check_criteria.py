# from criteria import criterias # for now - will need to be cleverer/cleaner/better 

from cel import evaluate

# print(criterias)




def construct_query(crt: dict):

    _filter = []
    for key in crt:

            

construct_query()

results = [
    {"drugsubclass": "Carbapenemase", "species":"Salmonella enterica"}, #this should be report
    {"drugsubclass": "Carbapenemases","gene":"blabla"}, # this will not report
    {"drugsubclass": "Carbapenemase (MBL)", "gene": "bla1", "species":"Salmonella enterica"}, # this should report
    {"drugsubclass": "Carbapenemase (MBL)", "gene": "bla1", "species":"Stenotrophomonas maltophilia"}, # this not should report
    {"drugsubclass": "Carbapenemase (MBL)", "gene": "bla2", "species":"Stenotrophomonas maltophilia"}, # this  should report

]