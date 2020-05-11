def get_data_from(data_file, start = 0, stop = -1):
    f = open(data_file, encoding = 'utf-8', mode = 'r')
    l_row = [[float(y) for y in x.strip().split()] for x in f.readlines()[start:stop]]
    l_column = []
    for i in range(len(l_row[0])):
        l_column.append([x[i] for x in l_row])
    return l_column

def get_colors_with(number):
    steps = int((number + 2) / 3)
    z = []
    z1 = double_color(steps)
    z = z + [(i[0], i[1], 0) for i in z1]
    z = z + [(0, i[0], i[1]) for i in z1]
    z = z + [(i[1], 0, i[0]) for i in z1]
    return z

def double_color(steps):
    z = []
    for i in range(steps):
        z.append((1 if i/steps < 0.5 else 1 - 2*(i/steps - 0.5), 
                  1 if i/steps > 0.5 else 2*i/steps))
    return z