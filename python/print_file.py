from input_file import *
import datetime
now = datetime.datetime.now()


## Functions which prints data to file

def print_fozard(filename, delta_time):
    fi = open(filename, "w")


    print("MODEL USING FOZARD INPUT DATA",file=fi)
    print(now.strftime("%Y-%m-%d %H:%M"),file=fi)  
    print("Time used(min): %.1f" % delta_time, file=fi)

    print("\nINPUT",file=fi)
    print("%i timesteps of %.2f seconds to calculate %.1f hours" % (N, dt*60, final_time/60) ,file=fi)
    print("Bulk concentration: %.2f" % conc_bulk,file=fi)
    print("Box size: %i, %i, %i" % (Lx/vl, Ly/vl, Lz/vl) ,file=fi)
    print("Side length: %i" % vl,file=fi)
    print("Positive feedback Kq: %i" % Kq,file=fi)
    fi.close

def print_output(filename, biofilm):
    dat = open(filename, "w")
    # xyz, cs, cqsm, cqsi
    print("x y z Np up down EPS subst qsm",file=dat)
    for v in biofilm.vortex_arr:
       print(v.get_pos(), v.get_num_particles(), v.cell_up_arr, v.cell_down_arr, v.eps_amount, v.conc_subst, v.conc_qsm,file=dat)

    dat.close()


