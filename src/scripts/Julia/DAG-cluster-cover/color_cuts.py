import re
import pydot
import glob

cuts_lst = glob.glob("*_cuts_*.txt") 
# j = 0
print(cuts_lst)
for i in cuts_lst:
    pattern_name = r"((wt\d{3}_\d{3}_\d_)cuts_(\d{2,}))\.txt"
    match_name = re.search(pattern_name,i)     
    cuts = open(i, 'r')

    j = 0

    for it in cuts.readlines():
        file_path = match_name.group(2) + "graph.dot"
        pattern_cuts = r"([-+])?\sx(\d+)"
        print(file_path)
        m_cuts = re.findall(pattern_cuts, it)
        i = 0
        file_handler = open(file_path, 'r')
        name_file = match_name.group(1) + "_represent_%03d" % (j)
        print(name_file)
        file_new = open(name_file + ".dot", 'w')
        for line in file_handler:
            color = "red" if m_cuts[i][0] == '-' else "green"
            a = int(m_cuts[i][1])
            b = str(a)
            pattern = r'(\d+) ->\s(\d+)\[label\s=\s+\"(%s)\"(\s?.+)\]' % (b)
            m = re.search(pattern, line)
            if m:
                file_new.write("%s -> %s [label=\"%s\" %s color=%s];\n" %
                               (m.group(1), m.group(2), m.group(3), m.group(4), color))
                i = i + 1
                if i == len(m_cuts):
                    i = 0
            else:
                file_new.write(line)

        j = j + 1

        file_handler.close()
        file_new.close()
        (graph,) = pydot.graph_from_dot_file(name_file + ".dot")
        graph.write_pdf(name_file + ".pdf")
    
    cuts.close()
