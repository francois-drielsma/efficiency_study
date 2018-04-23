import os
import shutil

class MergeCutsSummaryTex(object):
    def __init__(self):
        self.summary_list = []
        self.headings = []
        self.data = []
        self.caption = ""

    def append_summary(self, config, dir_selection):
        for i, a_dir in enumerate(config.conglomerate_dir):
            if i not in dir_selection:
                continue
            self.summary_list.append(a_dir+config.cuts_tex)
    
    def parse_one_line(self, line):
        words_tmp = line.split("&")
        words = []
        for word in words_tmp:
            #word = word.replace(" ", "")
            word = word.replace("//", "")
            word = word.replace("\n", "")
            words.append(word)
        return words

    def update_headings_list(self, line_number, words):
        if line_number >= len(self.headings):
            self.headings.append(words[0])
        elif self.headings[line_number] != words[0]:
            raise ValueError("could not match heading "+str(self.headings[line_number])+" to input "+str(words[0]))

    def update_data(self, file_number, line_number, words):
        row = line_number
        column = file_number

        if row == len(self.data):
            self.data.append([])
        if len(words) < 2:
            self.data[row].append(None)
        else:
            self.data[row].append(words[1])

    def print_data(self, folder, file_name):
        try:
            shutil.rmtree(folder)
        except OSError:
            pass
        os.makedirs(folder)
        fout = open(os.path.join(folder, file_name), "w")
        n_rows = len(self.data[0])
        head_matter = """
\\documentclass{letter}
\\usepackage{pdflscape}
\\begin{document}
"""
        tail_matter = """
\\end{document}
"""

        print >> fout, """
\\newcommand{\splitcell}[2][c]{%
    \\begin{tabular}[#1]{@{}c@{}}#2\end{tabular}}

\\begin{landscape}
\\begin{table}
\\centering
\\caption{"""+self.caption+"""}
\\begin{tabular}[pos]{l|""",
        for i in range(len(self.data[0])):
            fout.write("c")
        print >> fout, "}"
        for row in range(len(self.data)):
            #if row == 0:
            #    continue
            a_head = self.headings[row].replace("_", " ")
            print >> fout, a_head,
            for column in range(len(self.data[row])):
                if self.data[row][column] == None:
                    continue
                item = self.data[row][column]
                if row == 0:
                    item = "\\splitcell{"+item.replace(" ", "\\\\")+"}"
                print >> fout, "&", item,
            if self.headings[row] != "\hline":
                print >> fout, "\\\\",
            print >> fout
        print >> fout, """
\\end{tabular}
\\end{table}
\\end{landscape}
"""


    def latex(self, folder, file_name):
        return
        here = os.getcwd()
        os.chdir(folder)
        subprocess.check_output(["pdflatex", file_name])
        os.chdir(here)

    def merge_summaries(self, folder, file_name):
        for file_number, summary in enumerate(self.summary_list):
            print "Merging", summary
            try:
                fin = open(summary)
            except IOError:
                print "Failed to open", summary
                continue
            for line_number, line in enumerate(fin.readlines()):
                words = self.parse_one_line(line)
                if len(words) == 0:
                    continue
                self.update_headings_list(line_number, words)
                self.update_data(file_number, line_number, words)
        self.print_data(folder, file_name)
        self.latex(folder, file_name)
