def fasta_to_dictionaries():
    fast = ''
    while fast=='':
        try:
            link = input('\nEnter link to the fasta file of aligned sequences: ')
            #link = '/home/myshco/Завантажене/LS2.fasta'
            with open(link, 'r', encoding = 'utf8') as fast:
                falt = fast.read().strip('>').split('>')
        except:
            print("Couldn't read the file. Try again")
    num_of_data = len(falt)
    if num_of_data==1:
        print('Only one organism was detected. Unable to build a tree.')
        exit()
    print('Number of organisms: '+str(num_of_data))
    names={}
    seq={}
    for i in range(num_of_data):
        nl=int(falt[i].find('\n'))
        names[i]='>'+str(falt[i][0:nl])
        seq[i]=''.join(falt[i][nl+1:].split('\n'))
    return names, seq

def newick(link, nodes):
    l = len(nodes)
    line = ''
    line = new_node(nodes, l-1, line) + ';'
    file = open(link, 'w', encoding = 'utf8')
    file.write(line)
    file.close()
        
        
def new_node(nodes, n, line):
    if type(nodes[n][2]) == str:
        line = line + "('" + ''.join(nodes[n][2][1:].split(',')) + "':" + str(nodes[n][1]) + ','
    else:
        line = line + '('
        line = new_node(nodes, nodes[n][2][1], line) + ':' + str(nodes[n][1]) + ','
    if type(nodes[n][4]) == str:
        line = line + "'" + ''.join(nodes[n][4][1:].split(',')) + "':" + str(nodes[n][3]) + ')'
    else:
        line = new_node(nodes, nodes[n][4][1], line) + ':' + str(nodes[n][3]) + ')'
    return line



def d_to_matrix(names, seq):
    the_matrix=[]
    num = len(seq)
    for i in range(num):
        line = []
        line.append(names[i])
        a=seq[i]
        for j in range(num):
            b=seq[j]
            line.append(distance(a, b))
        the_matrix.append(line)
        
    return the_matrix

def distance(a, b):
    mismatches=0
    for f in range(len(a)):
        try:
            if not (a[f]==b[f]):
                mismatches+=1
        except:
            print('Sequences should be aligned')
            exit()
    return (mismatches)

def q_matrix(m):
    the_matrix=[]
    num = len(m)
    for i in range(num):
        line = []
        line.append(m[i][0])
        for j in range(num):
            em = 0
            en = 0
            for t in range(num):
                em += m[i][t+1]
            for y in range(num):
                en += m[j][y+1]
            if j==i:
                line.append(0.0)
            else:
                line.append((num-2) * m[i][j+1] - en - em)
        the_matrix.append(line)
    return the_matrix

def neighbour_joining(m, seq):
    nodes = []
    next_m = m
    l = len(m)
    for iteration in range(l-2):
        # Finding min
        q_m = q_matrix(next_m)
        minimum = 0
        for row in range(l-iteration):
            a_min = min(q_m[row][1:])
            for coli in range(l-iteration):
                if q_m[row][coli+1]==a_min:
                    column = coli
            if a_min < minimum:
                minimum = a_min
                row_of_min = row
                col_of_min = column
        # Calculating distance of pair members to node
        r_c_sum = 0
        for r in range(l-iteration):
            r_c_sum += (next_m[r][row_of_min+1] - next_m[r][col_of_min+1])
        row_d = (next_m[row_of_min][col_of_min+1] / 2) + r_c_sum / (2 * (l - iteration - 2))
        col_d = (next_m[row_of_min][col_of_min+1] / 2) - r_c_sum / (2 * (l - iteration - 2))
        nodes.append([iteration, row_d, next_m[row_of_min][0], col_d, next_m[col_of_min][0]])
        # Appending averaged distances
        try:
            a = int(next_m[row_of_min][0][0])
        except:
            a = 1
        try:
            b = int(next_m[col_of_min][0][0])
        except:
            b = 1
        new_row = [[a+b, iteration]] #name of clusters consists of number of leafs and № of node minus one (iteration)
        for row in range(l-iteration):
            new_distance = (next_m[row][row_of_min+1] + next_m[row][col_of_min+1] - next_m[row_of_min][col_of_min+1]) / 2
            next_m[row].append(new_distance)
            new_row.append(new_distance)
        new_row.append(0.0)
        next_m.append(new_row)
        # Deleting clustered elements
        next_m.pop(col_of_min)
        next_m.pop(row_of_min)
        for row in range(l-iteration-1):
            next_m[row].pop(col_of_min+1)
            next_m[row].pop(row_of_min+1)
    nodes.append([l-2, next_m[0][2], next_m[0][0], next_m[1][1], next_m[1][0]]) 
    return nodes


def mini(rist):
    minim = max(rist)
    for r in range(len(rist)):
        dat = rist[r]
        if (dat <= minim) and not (dat==0.0):
            minim = dat
            col = r
    return minim, col


def upgma(m, u_w):
    next_m = m
    l = len(m)
    nodes = []
    for iteration in range(l-1):
        # Finding min
        minimum = next_m[0][2]
        for row in range(l-iteration):
            a_min, colomn = mini(next_m[row][1:])
            if a_min <= minimum:
                minimum = a_min
                row_of_min = row
                col_of_min = colomn
        nodes.append([iteration, minimum/2, next_m[row_of_min][0], minimum/2, next_m[col_of_min][0]])
        # Appending averaged distances
        try:
            a = int(next_m[row_of_min][0][0])
            
        except:
            a = 1
        try:
            b = int(next_m[col_of_min][0][0])
        except:
            b = 1
        new_row = [[a+b, iteration]] #name of clusters consists of number of leafs and № of node minus one (iteration)
        for row in range(l-iteration):
            if u_w==True:
                new_distance = (next_m[row][row_of_min+1] * a + next_m[row][col_of_min+1] * b) / (a + b)
            else:
                new_distance = (next_m[row][row_of_min+1] + next_m[row][col_of_min+1]) / 2
            next_m[row].append(new_distance)
            new_row.append(new_distance)
        new_row.append(0.0)
        next_m.append(new_row)
        # Deleting clustered elements
        next_m.pop(row_of_min) #####################
        next_m.pop(col_of_min) #Do not swop col/row!
        for row in range(l-iteration-1):
            next_m[row].pop(row_of_min+1) #####################
            next_m[row].pop(col_of_min+1) #Do not swop col/row!
    li = len(nodes)
    for i in range(li):
        if type(nodes[li-i-1][2])==list:
            nodes[li-i-1][1] -= nodes[nodes[li-i-1][2][1]][1]
            
        if type(nodes[li-i-1][4])==list:
            nodes[li-i-1][3] -= nodes[nodes[li-i-1][4][1]][1]
    
    return nodes



def ordering(nodes, n, dots, e):
    a = nodes[n][2]
    if type(a)==list:
        dots, e = ordering(nodes, a[1], dots, e)
    else: #writing leaf down
        dots.append([e, a])
        e += 1
    
    dots.append([e, n])
    e += 1
    
    a = nodes[n][4]
    if type(a)==list:
        dots, e = ordering(nodes, a[1], dots, e)
    else: #writing leaf down
        dots.append([e, a])
        e += 1
    return dots, e

def fill_columns(nodes, af, n, bf, dots, tree):
    a = nodes[n][2]
    if type(a)==list:
        for s, name in dots:
            if a[1]==name:
                arrow = s
    else:
        for s, name in dots:
            if a==name:
                arrow = s
    
    b = nodes[n][4]
    if type(b)==list:
        for s, name in dots:
            if b[1]==name:
                brrow = s
    else:
        for s, name in dots:
            if b==name:
                brrow = s

    for s, name in dots:
            if n==name:
                nrrow = s
        
    for rrow in range(af, arrow):
        tree[rrow].append(' ')
    tree[arrow].append('┌')
    for rrow in range(arrow+1, nrrow):
        tree[rrow].append('│')
    tree[nrrow].append('┤')
    for rrow in range(nrrow+1, brrow):
        tree[rrow].append('│')
    tree[brrow].append('└')
    for rrow in range(brrow+1, bf):
        tree[rrow].append(' ')
    for br in range(int(nodes[n][1]*resolution)):
        for rrow in range(af, arrow):
            tree[rrow].append(' ')
        tree[arrow].append('─')
        for rrow in range(arrow+1, nrrow):
            tree[rrow].append(' ')
            
    for br in range(int(nodes[n][3]*resolution)):
        for rrow in range(nrrow+1, brrow):
            tree[rrow].append(' ')
        tree[brrow].append('─')
        for rrow in range(brrow+1, bf):
            tree[rrow].append(' ')
    
    if type(a)==list:
        tree = fill_columns(nodes, af, a[1], nrrow, dots, tree)
    else:
        tree[arrow].append(dots[arrow][1])
    if type(b)==list:
        tree = fill_columns(nodes, nrrow, b[1], bf, dots, tree)
    else:
        tree[brrow].append(dots[brrow][1])
        
    return tree
    
    
def drawing_a_tree(nodes, resolution):
    l=len(nodes)
    # Finding order
    dots = []
    e = 0
    dots, e = ordering(nodes, (l-1), dots, e)
    # Filling colomns
    tree = []
    for n in dots:
        tree.append([' '])
    af = 0
    bf = len(tree)
    tree = fill_columns(nodes, af, (l-1), bf, dots, tree)
    the_tree = []
    for i in range(len(tree)):
        the_tree.append(''.join(tree[i]))
    str_tree = '\n'.join(the_tree)
    print(str_tree)
#│┌└ ┤─



print('''Welcome to the tree calculator.
Build UPGMA, WPGMA or neighbour joining trees from aligned nucleotides/aminoacids sequence.''')
d_names, d_seq = fasta_to_dictionaries()
matrix = d_to_matrix(d_names, d_seq)

nodes_of_tree = []
while nodes_of_tree==[]:
    algorithm = input('Choose algorithm (UPGMA/WPGMA/NJ): ')
    if algorithm.upper()=='UPGMA':
        print('Calculating UPGMA tree...')
        nodes_of_tree = upgma(matrix, True)
        print('Done')
    elif algorithm.upper()=='WPGMA':
        print('Calculating WPGMA tree...')
        nodes_of_tree = upgma(matrix, False)
        print('Done')
    elif algorithm.upper()=='NJ' or algorithm.upper()=='NEIGHBOUR JOINING':
        print('Calculating neighbour joining tree...')
        nodes_of_tree = neighbour_joining(matrix, d_seq)
        print('Done')
    else:
        print('Please, choose one of the methods')

do_newick = input('\nWould you like to save your tree as Newick file? (Y/n): ')
if do_newick=='Y':
    saving_link = input('Choose file name and write it with destination folder (like link you entered to fasta file): ')
    newick(saving_link, nodes_of_tree)
    print('Your file have been saved to the folder. You can open it using https://www.chiplot.online/normalTree.html or software you like.\nnormalTree can take some time to load')

do_drawing = input('\nWould you like to print the tree to the shell? (Y/n): ')
if do_drawing=='Y':
    resolution = 1
    drawing_a_tree(nodes_of_tree, resolution)
    print('''\n
        Try zooming out (Ctrl, -) or in (Ctrl, Shift, +) if you cannot see a nice tree.
        You can try changing the font of your shell if you want.
        If you used UPGMA or WPGMA, short branch lengths can be equal due to minimal resolution of one character.
        If there are problems, try opening Newick file in more powerfull program. You may use https://www.chiplot.online/normalTree.html in your web browser.''')
        
print('\nThank you for using our calculator!')