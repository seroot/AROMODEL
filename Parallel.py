#! /usr/bin/python

import Configure

def Write_Submit_Script(Num_Atoms, Sim_Name, Name, GPU = False):
    NProcs = int(Num_Atoms/1000.)
    NProcs = 24
    if NProcs > 48:
        GPU = True
    # Set up submit script
    if not GPU:
        if NProcs <= 24:
            sub_temp = Configure.Template_Path + "sub_Lammps"
            submit = "sub_%s" % Name
            with open(sub_temp) as f:
                template = f.read()
            s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % Name, NProcs = NProcs, Nodes=1, tpn = NProcs)
            with open(submit,'w') as f:
                f.write(s)
        if NProcs > 24 and NProcs < 48:
            sub_temp = Configure.Template_Path + "sub_Lammps"
            submit = "sub_%s" % Name
            tpn = NProcs/2
            NProcs = tpn*2
            with open(sub_temp) as f:
                template = f.read()
            s = template.format(Sim_Name= Sim_Name, path = Configure.Comet_Path % Name, NProcs = NProcs, Nodes=2, tpn = tpn)
            with open(submit,'w') as f:
                f.write(s)

    elif GPU:
        sub_temp = Configure.Template_Path + "GPU_Sub"
        submit = "GPU_sub_%s" % Name
        with open(sub_temp) as f:
            template = f.read()
        s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % Name)
        with open(submit,'w') as f:
            f.write(s)

    return submit