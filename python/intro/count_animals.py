def count_moose(file):
    '''
    This function counts the number of moose
    '''
    newfile = open(file, 'r') 
    filelines = newfile.readlines()

    moose_list = []
    for line in filelines:
        d, t, a, n = line.split(' ')
        if a == 'Moose':
            num = int(n)
            if num > 10:
                moose_list.append(line)

    return moose_list
