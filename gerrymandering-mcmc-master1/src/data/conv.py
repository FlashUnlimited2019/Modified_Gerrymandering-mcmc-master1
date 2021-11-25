import json
a_file = open(
    r'src\data\WA_final.json', "r")
json_object = json.load(a_file)
a_file.close()
print(json_object['features'][0]['properties']['POP100'])
