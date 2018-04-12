#!/d/home/apps/build64/bin/python

"""
 .. This GUI is meant to intialize the NTA Data Processing Packages
    Module Author: Hussein Al Ghoul hussein.al-ghoul@epa.gov 
"""
import Tkinter as tk
import Tkconstants, tkFileDialog
import tkMessageBox
import numpy, os
import functions_Universal_v1 as fn
import batch_search_v1 as bs
import Toxpi_v1 as tp
from threading import Thread
#from graphicHelp import help
import time

DIRECTORY = None
ERROR = 1
OUTPUT = ""


class Control(object):

    df=[0,0]
    dft = None
    Fname = ['','']
    Sname = ['','']
    Rname = ['','']
    name = ['','']
    Files = ['','']
    def __int__(self, study = None, user = None):

        self.studyL = study
        self.userL = user


    def study(self):
        st = ent0.get()
        self.studyL = st
    
    def user(self):
        us = ent1.get()
        self.userL = us


    def run(self):
    
        self.study()
        self.user()

    def load(self,con):

        self.studyL = con[0]
        self.userL = con[1]

        ent0.insert("end",con[0])
        ent1.insert("end",con[1])


    def check_file_exists(arg,f,index):
        answer = ['','']
        file = ['','']
        file[index]=f
        FFPath = ['','']
        message = file[index] + " Exists. Would You Like to Replace It?"
        FFPath[index]=os.getcwd()+"/"+DIRECTORY+"/"+ file[index]
        timestr = time.strftime('%Y%m%d')
        if os.path.exists(FFPath[index]):
            file[index] = file[index].rsplit('.',1)[0] + "_" + timestr + ".csv"

        print file[index]
        return file[index]


    def initialize(self):
        self.run()
        global DIRECTORY
        con = [self.studyL,self.userL]
        numpy.save(os.getcwd()+"/control_list.npy",con)
        if ERROR == 0:
            #frame2.grid(row=0,column=2,sticky="N"+"S"+"E"+"W")
            frame2.pack(side='left',expand=1,fill='both')
            btnrd.pack(side="left")
            btndd.pack(side="left")
            ckbxrd3.pack(side='left')
            btnstat.pack(side="left")
            ckbxrd4.pack(side='left')
            btnct.pack(side="left")
            ckbxrd31.pack(side='left')
            btnbk1.pack(side="left")
            btnxt1.pack(side="left")
            root.update()
            if ERROR == 1:
                tkMessageBox.showinfo("ERROR","You Need to Select a Data File")
                DIRECTORY = con[0] + "_" + con[1]
                if not os.path.exists(DIRECTORY):
                    os.makedirs(DIRECTORY)
        return DIRECTORY


    def openfile(self):
        global ERROR
        con = [self.studyL,self.userL]
        ent2.delete(0,tk.END)
        numpy.save(os.getcwd()+"/control_list.npy",con)
        root.filename = tkFileDialog.askopenfilenames(initialdir = os.getcwd(),title = "Select file",filetypes = (("csv files","*.csv"),("tsv files","*.tsv"),("all files","*.*")))
        fname = ['','']
        for i in range(2):
            fname[i] = root.tk.splitlist(root.filename)[i].rsplit('/',1)[-1]
        fname_list = fname[0] +";"+fname[1]
        contra.Files = root.tk.splitlist(root.filename)            
        ent2.insert("end",fname_list)
        ERROR = 0
        return contra.Files
    #return root.filename

    def next1(self):
        frame3.pack(expand=1,fill='both')
        #frame3.grid(row=0,column=3,sticky="N"+"S"+"E"+"W")
        btncf.pack(side="left")
        ckbxrd32.pack(side='left')
        btncflags.pack(side="left")
        ckbxrd33.pack(side='left')
        btncd.pack(side="left")
        ckbxrd34.pack(side='left')
        btntp.pack(side="left")
        ckbxrd35.pack(side='left')
        btnbk2.pack(side="left")     
        btnxt2.pack(side="left")
        root.update()


    def next2(self):
        frame4.grid(row=0,column=3,sticky="N"+"S"+"E"+"W")
        btnAccNI.pack(side="top")
        btnDLAccNI.pack(side="top")
        root.update()


    def back1(self):
        frame2.pack_forget()
        root.update()


    def back2(self):
        frame3.pack_forget()
        root.update()


    def Read_Data(arg,File,index):
        global OUTPUT
        btnrd.configure(bg='light yellow')
        contra.df[index] = fn.read_data(File,index)
        OUTPUT = "Data was converted to Dataframe \n"
        if index == 1:
            btnrd.configure(bg='pale green')
        T.insert(tk.END, OUTPUT)
        return contra.df[index]


    def Parse_Data(arg,index):
        global OUTPUT
        print contra.df[index]
        contra.df[index] = fn.parse_headers(contra.df[index],index)
        OUTPUT= "Data File was parsed! \n" 
        if index == 1:
            btnpd.configure(bg='pale green')
        T.insert(tk.END, OUTPUT)
        return contra.df[index]      


    def Statistics(arg,File,index):
        global OUTPUT
        btnstat.configure(bg='light goldenrod yellow')
        contra.Fname[index] = File.rsplit('/',1)[-1]
        contra.Sname[index] = contra.Fname[index].rsplit('.',1)[0] + "_Statistics.csv"
        contra.df[index] = fn.statistics(contra.df[index],index)
        if checkCmd2.get():
            OUTPUT = "Will Save the Statistics File \n"
            contra.name[index]=contra.check_file_exists(contra.Sname[index],index)
            contra.df[index].to_csv(os.getcwd()+"/"+DIRECTORY+"/"+contra.name[index], index=False)
        else:
            OUTPUT = "Not Saving the Statistics File \n"
        T.insert(tk.END, OUTPUT)
        if index == 1:
            btnstat.configure(bg='pale green')
        OUTPUT = "Statistical Parsing is done! \n"
        return contra.df[index]      


    def Check_Tracers(arg,File,index):
        global OUTPUT
        btnct.configure(bg='light yellow')
        df_tracers = [None,None]
        contra.Fname[index] = File.rsplit('/',1)[-1]
        contra.Sname[index] = contra.Fname[index].rsplit('.',1)[0] + "_Tracers.csv"
        df_tracers[index] = fn.check_feature_tracers(contra.df[index],index,10,0.5)
        if checkCmd3.get():
            OUTPUT = "Will Save the Tracers File \n"
            contra.name[index]=contra.check_file_exists(contra.Sname[index],index)
            df_tracers[index].to_csv(os.getcwd()+"/"+DIRECTORY+"/"+contra.name[index], index=False)
        else:
            OUTPUT = "Not Saving the Tracers File \n"
        T.insert(tk.END, OUTPUT)
        if index == 1:
            btnct.configure(bg='pale green')
        OUTPUT = "Tracers check is done! \n"
        T.insert(tk.END, OUTPUT)
        return contra.df[index]      


    def Clean_Features(arg,File,index):
        global OUPUT
        btncf.configure(bg='light yellow')
        contra.Fname[index] = File.rsplit('/',1)[-1]
        contra.Sname[index] = contra.Fname[index].rsplit('.',1)[0] + "_Clean.csv"
        contra.df[index] = fn.clean_features(contra.df[index],index)
        if checkCmd4.get():
            OUTPUT = "Will Save the Cleaned File \n"
            contra.name[index]=contra.check_file_exists(contra.Sname[index],index)
            contra.df[index].to_csv(os.getcwd()+"/"+DIRECTORY+"/"+contra.name[index], index=False)
        else:
            OUTPUT = "Not Saving the Cleaned File \n"
        T.insert(tk.END, OUTPUT)
        contra.df[index] = fn.Blank_Subtract(contra.df[index],index)
        contra.df[index].to_csv(os.getcwd()+"/"+DIRECTORY+"/"+contra.Fname[index].rsplit('.',1)[0]+"_reduced.csv", index=False)
        if index == 1:
            btncf.configure(bg='pale green')
        OUTPUT = "Features were cleaned! \n"
        T.insert(tk.END, OUTPUT)
        return contra.df[index]



    def Create_Flags(arg,File,index):
        global OUPUT
        btncflags.configure(bg='light yellow')
        contra.Fname[index] = File.rsplit('/',1)[-1]
        contra.Sname[index] = contra.Fname[index].rsplit('.',1)[0] + "_Flags.csv"
        contra.Rname[index] = contra.Fname[index].rsplit('.',1)[0] + "_Reduced.csv"    
        contra.df[index] = fn.flags(contra.df[index])
        if checkCmd5.get():
            OUTPUT = "Will Save the Flags File \n"
            contra.name[index]=contra.check_file_exists(contra.Sname[index],index)
            contra.df[index].to_csv(os.getcwd()+"/"+DIRECTORY+"/"+contra.name[index], index=False)
        #fn.reduce(contra.df[index],index).to_csv(os.getcwd()+"/"+DIRECTORY+"/"+contra.Rname[index], index=False)
        else:
            OUTPUT = "Not Saving the Flags File \n"
        T.insert(tk.END, OUTPUT)
        if index == 1:
            btncflags.configure(bg='pale green')
        OUTPUT = "Flags were created! \n"
        T.insert(tk.END, OUTPUT)
        return contra.df[index]


    def Drop_Duplicates(arg,File,index):
        global OUPUT
        btndd.configure(bg='light yellow')
        contra.Fname[index] = File.rsplit('/',1)[-1]
        contra.Sname[index] = contra.Fname[index].rsplit('.',1)[0] + "_AfterDuplicates.csv"
        contra.df[index] = fn.duplicates(contra.df[index],index)
        if checkCmd1.get():
            OUTPUT = "Will Save the Duplicates File \n"
            contra.name[index]=contra.check_file_exists(contra.Sname[index],index)
            contra.df[index].to_csv(os.getcwd()+"/"+DIRECTORY+"/"+contra.name[index], index=False)
        else:
            OUTPUT = "Not Saving the Duplicates File \n"
        if index == 1:
            btndd.configure(bg='pale green')
        T.insert(tk.END, OUTPUT)
        OUTPUT = "Duplicates were Removed! \n"
        T.insert(tk.END, OUTPUT)
        return contra.df[index]      


    def Combine_Modes(self):
        global OUPUT
        btncd.configure(bg='light yellow')
        contra.Fname[0] = contra.Files[0].rsplit('/',1)[-1]
        contra.Sname[0] = contra.Fname[0].rsplit('.',1)[0] + "_Combined.csv"
        contra.dft = fn.combine(contra.df[0],contra.df[1])
        if checkCmd6.get():
            OUTPUT = "Will Save the Combined File \n"
            contra.name[0]=contra.check_file_exists(contra.Sname[0],0)
            contra.dft.to_csv(os.getcwd()+"/"+DIRECTORY+"/"+contra.name[0], index=False)
        else:
            OUTPUT = "Not Saving the Combined File \n"
        btncd.configure(bg='pale green')
        directory = os.getcwd()+"/"+DIRECTORY
        fn.MPP_Ready(contra.dft,directory,contra.Fname[0].rsplit('.',1)[0])
        T.insert(tk.END, OUTPUT)
        OUTPUT = "DModes were Combined! \n"
        T.insert(tk.END, OUTPUT)
        return contra.dft 


    def Tox_Pi(self):
        global OUPUT
        btntp.configure(bg='light yellow')
        contra.Batch_Search()
        directory = os.getcwd()+"/"+DIRECTORY
        contra.Fname[0] = contra.Files[0].rsplit('/',1)[-1]
        contra.Sname[0] = contra.Fname[0].rsplit('.',1)[0] + "_toxpi.csv"
        print "Finished the Selenium part"
        time.sleep(10)
    #dashboard_file[0] = [filename for filename in os.listdir(os.getcwd()+"/"+DIRECTORY) if filename.startswith('ChemistryDashboard-AdvancedSearch')]
        dashboard_file = contra.Download_Finished()
        print "This is the dashboard_file: " + dashboard_file
        contra.dft = tp.process_toxpi(contra.dft,directory, dashboard_file)
        contra.dft = tp.calculate_toxpi(contra.dft,directory)
        if checkCmd7.get():
            OUTPUT = "Will Save the ToxPi File \n"
            contra.name[0]=contra.check_file_exists(contra.Sname[0],0)
            contra.dft.to_csv(os.getcwd()+"/"+DIRECTORY+"/"+contra.name[0], index=False)
        else:
            OUTPUT = "Not Saving the toxpi File \n"
        btntp.configure(bg='pale green')
        T.insert(tk.END, OUTPUT)
        OUTPUT = "Toxpi Data was Created! \n"
        T.insert(tk.END, OUTPUT)



    def plot_window(self):
        plotw = tk.Toplevel()
        plotw.title("ToxPi")
        btnplot=tk.Button(root,text="Plot" , height=2 , width = 15,command = lambda: contra.plot())

    def Batch_Search(self):
        directory = os.getcwd()+"/"+DIRECTORY
        compounds = []
        compounds = fn.formulas(contra.dft)
        bs.batch_search(compounds[:100],directory)
    #bs.batch_search(compounds,directory)

    def Download_Finished(self):
        directory = os.getcwd()+"/"+DIRECTORY
        finished = False
        file = None
        while not finished:
            for filename in os.listdir(directory):
                if not filename.startswith('ChemistryDashboard-Batch-Search'):
                    time.sleep(2)
                if filename.startswith('ChemistryDashboard-Batch-Search'):
                    print "in loop filename: "+filename
                    file = filename
                    finished = True
        print "This is what was downloaded: " + file
        return file
    

# Multiprocessing Functions

    def RD_MP(self):
        threads = []    
        for index, f in enumerate(contra.Files):
            thread = Thread(target=contra.Read_Data,args=(f,index))
            threads.append(thread)
        #thread.start()
        for thread in threads:
            thread.start()



    def DD_MP(self):
        threads = []
        for index, f in enumerate(contra.Files):
            thread = Thread(target=contra.Drop_Duplicates,args=(f,index))
            threads.append(thread)
        for thread in threads:
            thread.start()


    def S_MP(self):
        threads = []
        for index, f in enumerate(contra.Files):
            thread = Thread(target=contra.Statistics,args=(f,index))
            threads.append(thread)
        for thread in threads:
            thread.start()


    def CF_MP(self):
        threads = []
        for index, f in enumerate(contra.Files):
            thread = Thread(target=contra.Clean_Features,args=(f,index))
            threads.append(thread)
        for thread in threads:
            thread.start()


    def F_MP(self):
        threads = []
        for index, f in enumerate(contra.Files):
            thread = Thread(target=contra.Create_Flags,args=(f,index))
            threads.append(thread)
        for thread in threads:
            thread.start()


    def CT_MP(self):
        threads = []
        for index, f in enumerate(contra.Files):
            thread = Thread(target=contra.Check_Tracers,args=(f,index))
            threads.append(thread)
        for thread in threads:
            thread.start()

    def BS_MP(self):
        threads = []
        for index, f in enumerate(contra.Files):
            thread = Thread(target=contra.Blank_Subtract,args=(f,index))
            threads.append(thread)
        for thread in threads:
            thread.start()


lst = ['STUDY','USER','DATA']



# Main GUI window
contra = Control()
root = tk.Tk()
root.wm_title("NTA INITIALIZER")

frame=tk.LabelFrame(root, text = 'Details' ,height = 600, width = 600)
frame2=tk.LabelFrame(root, text = 'Step One' ,height = 600, width = 600)
frame3=tk.LabelFrame(root, text = 'Step Two' ,height = 600, width = 600)
frame4=tk.LabelFrame(root, text = 'Step Three' ,height = 600, width = 600)

frame.pack(side='left',padx=10,expand=1,fill='both')


frame00=tk.Frame(frame)
frame00.pack(expand=1,fill='both')
l0=tk.Label(frame00, text=lst[0])
l0.pack(padx=5,pady=15,expand=1,side="left",fill='both')
ent0=tk.Entry(frame00, width=30, font=55)
ent0.pack(padx=5,pady=15,expand=1,side="left",fill='both')


frame01=tk.Frame(frame)
frame01.pack(expand=1,fill='both')
l1=tk.Label(frame01, text=lst[1])
l1.pack(padx=10,pady=15,expand=1,side="left",fill='both')
ent1=tk.Entry(frame01, width=30, font=55)
ent1.pack(padx=5,pady=15,expand=1,side="left",fill='both')


frame02=tk.Frame(frame)
frame02.pack(expand=1,fill='both')
l2=tk.Label(frame02, text=lst[2])
l2.pack(padx=5,pady=15,expand=1,side="left",fill='both')
ent2=tk.Entry(frame02, width=22, font=15)
ent2.pack(padx=5,pady=15,expand=1,side="left",fill='both')
btnfile=tk.Button(frame02, text="Open",  height=1 , width = 5,command = lambda: contra.openfile())
btnfile.pack(padx=5,pady=15,expand=0,side="left",fill='both')



frame1=tk.Frame(frame)
frame1.pack(expand=1,side="top")

btnSt=tk.Button(frame1, text="Initialize" , command= lambda: contra.initialize())
btnSt.pack(expand=1,side='top')

####first frame where binning is done#####
checkCmd1 = tk.IntVar()
checkCmd2 = tk.IntVar()
checkCmd3 = tk.IntVar()

frame21=tk.Frame(frame2)
frame21.pack(padx=80,expand=1,fill='both')
btnrd=tk.Button(frame21, text="Read Data" ,  height=2 , width = 15,command = lambda: contra.RD_MP())


frame22=tk.Frame(frame2)
frame22.pack(padx=80,expand=1,fill='both')
btndd=tk.Button(frame22,text="Drop Duplicates" , height=2 , width = 15,command = lambda: contra.DD_MP())
ckbxrd3=tk.Checkbutton(frame22,variable=checkCmd1 ,text="Save File")

frame23=tk.Frame(frame2)
frame23.pack(padx=80,expand=1,fill='both')
btnstat=tk.Button(frame23, text="Statistics" ,  height=2 , width = 15,command= lambda: contra.S_MP())
ckbxrd4=tk.Checkbutton(frame23,variable=checkCmd2 ,text="Save File")


frame24=tk.Frame(frame2)
frame24.pack(padx=80,expand=1,fill='both')
btnct=tk.Button(frame24, text="Tracers Check" ,  height=2 , width = 15,command= lambda: contra.CT_MP())
ckbxrd31=tk.Checkbutton(frame24,variable=checkCmd3 ,text="Save File")



btnxt1=tk.Button(frame2, text="Next",  height=2 , width = 15,command = lambda: contra.next1())
btnbk1=tk.Button(frame2, text="Back", height=2 , width = 15 ,command = lambda: contra.back1())



####second frame where Decay Amps are calculated####

checkCmd4 = tk.IntVar()
checkCmd5 = tk.IntVar()
checkCmd6 = tk.IntVar()
checkCmd7 = tk.IntVar()

frame31=tk.Frame(frame3)
frame31.pack(padx=80,expand=1,fill='both')
btncf=tk.Button(frame31,text="Clean Features" , height=2 , width = 15,command = lambda: contra.CF_MP())
ckbxrd32=tk.Checkbutton(frame31,variable=checkCmd4 ,text="Save File")

frame32=tk.Frame(frame3)
frame32.pack(padx=80,expand=1,fill='both')
btncflags=tk.Button(frame32, text="Create Flags", height=2, width =15, command = lambda: contra.F_MP())
ckbxrd33=tk.Checkbutton(frame32, variable=checkCmd5 ,text="Save File")

frame33=tk.Frame(frame3)
frame33.pack(padx=80,expand=1,fill='both')
btncd=tk.Button(frame33, text="Combine Modes" ,  height=2 , width = 15,command= lambda: contra.Combine_Modes())
ckbxrd34=tk.Checkbutton(frame33,variable=checkCmd6 ,text="Save File")

frame34=tk.Frame(frame3)
frame34.pack(padx=80,expand=1,fill='both')
btntp=tk.Button(frame34, text="ToxPi" ,  height=2 , width = 15,command= lambda: contra.Tox_Pi())
ckbxrd35=tk.Checkbutton(frame34,variable=checkCmd7 ,text="Save File")


btnxt2=tk.Button(frame3, text="Next",  height=2 , width = 15,command = lambda: contra.next2())
btnbk2=tk.Button(frame3, text="Back", height=2 , width = 15 ,command = lambda: contra.back2())


##### Toxpi plot window ######



# Output Message log Window
rootS = tk.Toplevel()
rootS.title("MESSAGE LOG")
S = tk.Scrollbar(rootS)
T = tk.Text(rootS , height=10, width=50)
S.pack(side='right', fill='both')
T.pack(side='left', fill='both')
S.config(command = T.yview)
T.config(yscrollcommand=S.set)
'''wm=300
hm=300
ws = rootS.winfo_screenwidth()
hs = rootS.winfo_screenheight()
xm = 0.45*ws
ym = 0.5*hs'''
#rootS.geometry('%dx%d+%d+%d' % (wm,hm,xm,ym))

hidden0 = tk.Toplevel()
hidden0.withdraw()

hidden1 = tk.Toplevel()
hidden1.withdraw()


if os.path.isfile(os.getcwd()+"/control_list.npy"):
    con = numpy.load(os.getcwd()+"/control_list.npy")
    contra.load(con)

#   print "Hi there"

root.mainloop()
#rootS.mainloop()

