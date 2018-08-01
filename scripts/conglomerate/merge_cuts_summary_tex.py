import glob
import os
import shutil
import utilities.cut_names

class MergeCutsSummaryTex(object):
    def __init__(self):
        # list of folders; each containing a list of summary documents
        self.summary_list_of_lists = []
        # builds a list of row headings
        self.headings = []
        # builds a list of row data
        self.data = []
        # table captions
        self.caption = ""
        self.vertical_splits = 8
        self.table_ref_pre = ""

    def append_summary(self, config, dir_selection):
        for i, a_dir in enumerate(config.conglomerate_dir):
            print "Adding...", a_dir+config.cuts_tex,
            if i not in dir_selection:
                continue
            self.summary_list_of_lists.append(glob.glob(a_dir+config.cuts_tex))
            print len(self.summary_list_of_lists), "files"
    
    def parse_one_line(self, line):
        words_tmp = line.split("&")
        words = []
        for word in words_tmp:
            #word = word.replace(" ", "")
            word = word.replace("//", "")
            word = word.replace("\n", "")
            words.append(word)
        return words

    def update_headings_list(self, file_number, line_number, words):
        if file_number >= len(self.headings):
            self.headings.append([])
        a_heading = words[0].rstrip(' ')
        if a_heading not in ["\\hline"]:
            a_heading = utilities.cut_names.cut_names[a_heading]
        if line_number >= len(self.headings[file_number]):
            self.headings[file_number].append(a_heading)
        elif self.headings[file_number][line_number] != a_heading:
            print "file", file_number, "line", line_number, "heading", a_heading
            for head in self.headings:
                print head
            raise ValueError("could not match heading "+str(self.headings[file_number][line_number])+" to input "+str(words[0]))

    def update_data(self, dir_number, file_number, line_number, words):
        row = line_number
        column = dir_number

        if file_number == len(self.data):
            self.data.append([])
        if row == len(self.data[file_number]):
            self.data[file_number].append([])
        if len(words) < 2:
            self.data[file_number][row].append(None)
        else:
            self.data[file_number][row].append(words[1])

    def print_data(self, file_number, folder, file_name):
        """
        Write merged data to disk
        - file_number: indexes the file from which the cuts data is taken
        - folder: is the directory to which the merged table is written
        - file_name: is the name to which the merged table is written; the bit 
                  before the .tex is used as the table reference \ref{tab:file_name}
        """
        fout = open(os.path.join(folder, file_name), "w")
        n_rows = len(self.data[file_number][0])
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
"""
        n_tables = (len(self.data[0][0])-1)/self.vertical_splits + 1
        for table_index in range(n_tables):
            start = table_index*self.vertical_splits
            end = min((table_index+1)*self.vertical_splits, len(self.data[0][0]))
            cell_format = "{l|"+"c"*(end-start)+"}\n"
            table_ref = file_name.split(".")[0]
            table_ref = self.table_ref_pre+table_ref+"_"+str(table_index)
            print >> fout, """
\\begin{landscape}
\\begin{table}
\\centering
\\caption{"""+self.caption[file_number][table_index]+"""\label{tab:"""+table_ref+"""}}
\\begin{tabular}[pos]"""+cell_format,
            for row in range(len(self.data[file_number])):
                #if row == 0:
                #    continue
                a_head = self.headings[file_number][row]
                print >> fout, a_head.ljust(utilities.cut_names.max_length),
                for column in range(start, end):
                    if self.data[file_number][row][column] == None:
                        continue
                    item = self.data[file_number][row][column]
                    if row == 0:
                        item = "\\splitcell{"+item.replace(" ", "\\\\")+"}"
                        item = item.replace("_", " ")
                    print >> fout, "&", str(item).rjust(8),
                if self.headings[file_number][row] != "\hline":
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

    def clean_dir(self, folder):
        try:
            shutil.rmtree(folder)
        except OSError:
            pass
        os.makedirs(folder)


    def merge_summaries(self, folder, file_prefix):
        for dir_number, summary_list in enumerate(self.summary_list_of_lists):
            for file_number, summary in enumerate(summary_list):
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
                    self.update_headings_list(file_number, line_number, words)
                    self.update_data(dir_number, file_number, line_number, words)
        self.clean_dir(folder)
        for file_number, summary in enumerate(self.summary_list_of_lists[0]):
            file_name = file_prefix+"_"+str(file_number)+".tex"
            self.print_data(file_number, folder, file_name)
            self.latex(folder, file_name)

    name_lookup = {
      
    }

