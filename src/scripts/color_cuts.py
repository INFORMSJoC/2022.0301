import re
import pydot

cuts = open("cuts_3/res_cuts_3_uniq.txt", 'r')

j = 0

for it in cuts.readlines():
    file_path = "graph.dot"
    pattern_cuts = r"([-+])?\sx(\d+)"
    m_cuts = re.findall(pattern_cuts, it)
    i = 0
    file_handler = open(file_path, 'r')
    name_file = "cut_represent_%d" % (j)
    file_new = open(name_file + ".dot", 'w')
    for line in file_handler:
        color = "red" if m_cuts[i][0] == '-' else "green"
        pattern = r'(\d+) ->\s(\d+)\[label\s=\s+\"(%s)\"(\s?.+)\]' % (m_cuts[i][1])
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
    graph.write_png(name_file + ".png")
