#Read input file
#   parameters:#       input file - .i input file
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
def read_input_file(input_file):
    inputs = {}
    with open(input_file) as f:
        for line in f:
            if '=' in line:
                inputs[line.split("=")[0].strip()] = line.split("=")[1].strip()
            else: pass
        for key,val in inputs.items():
            if is_number(val) == True and key != 'ObsIDs':
                inputs[key] = float(val)
            if key == 'ObsIDs':
                #Obtain individual obsids from list
                obsids = [inputs['ObsIDs'].split(',')[i].strip() for i in range(len(inputs['ObsIDs'].split(',')))]
                inputs['ObsIDs'] = obsids
    return inputs
