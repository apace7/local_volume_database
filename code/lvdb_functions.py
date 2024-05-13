import yaml

def get_notes(key):
    path="/Users/apace/Documents/local_volume_database/data_input/"
    with open(path+ key + '.yaml', 'r') as stream:
        try:
            stream_yaml = yaml.load(stream, Loader=yaml.Loader)
            if 'notes' in stream_yaml.keys():
                print("notes:\t", key)
#                 print(stream_yaml['notes'])
                for i in range(len(stream_yaml['notes'])):
                    print(stream_yaml['notes'][i])
            else:
                print("no notes", key)
        except:
            print("no key: ", key)