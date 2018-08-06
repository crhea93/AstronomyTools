#Read input file
#   parameters:
#       input file - .i input file
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
                inputs[line.split("=")[0].strip().lower()] = line.split("=")[1].strip()
            else: pass
        if len(inputs) != 11:
            print("Please recheck the input file since some parameter is missing...")
            print("Exiting program...")
            exit()
        else:
            print("Successfully read in input file")
            for key,val in inputs.items():
                if is_number(val) == True:
                    inputs[key] = float(val)
        return inputs
