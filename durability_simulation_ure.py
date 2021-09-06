#-------------------------------------------------------------------------------
# Name:        durability_simulation with unrecoverable errors
# Purpose:      Monte Carlo simualtion to estimate the probability of data loss in erasure coded systems
# Author:      TQ, quarktetra@gmail.com
# Created:     25/08/2021
# Requires:    Python3
#-------------------------------------------------------------------------------
# run this from the cmd window as: durability_simulation.py 20 2 1 1 20 50 1
# in IDE, use command line parameters: 20 2 1 1 20 50 1
# see https://tetraquark.netlify.app/post/raid_durability/   for details
import argparse
import math
import numpy as np
import random

def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n - 1)


def sifter(systems,p_shards,repair_time,failure_rate,hval):   # feed this with the surviving systems. It will return the instances of data losses, and the systems that might still lose data
    instances=0
    instances_wuer=0
    systems_to_keep=[]
    for sysid in range(0,len(systems)):
        failuretimes=systems[sysid]
        failuretimes=failuretimes[failuretimes<365.25]  # we are interested in the first year
        sys_lost_data=0
        sys_lost_data_wuer=0
        if failuretimes.size>p_shards-1: # for data loss you need to have more than p_shards failures, otherwise the system will not lose any data
            failuretimes=np.sort (failuretimes)   # time order the failures for each system. We will carefully count how many will accumulate over time
            toberecovered=np.array([failuretimes[0]])    # the first failure enters into the toberecovered bucket.
            for findex in range(1,failuretimes.size):
                this_time= failuretimes[findex]
                toberecovered=toberecovered[toberecovered>this_time-repair_time]   # keep the ones that have not been repaired yet. (this_time< t_failure+repair_time )
                toberecovered=np.append(toberecovered,[this_time])  # add the current failure to the toberecovered list.

                if len(toberecovered)>p_shards: # at any give point, if the list contains more items than the parity count, data is lost
                    sys_lost_data=1
                    sys_lost_data_wuer=1
                if len(toberecovered)>p_shards-1 and random.random()<=hval:   # pull a random number. If it is smaller than hval, then it is URE failures
                    sys_lost_data_wuer=1


            if sys_lost_data_wuer==0:
                failuretimes[0]=failuretimes[0]+np.random.exponential(failure_rate, 1) #the first failed drive will be replaced. Create a new failure time for it. (add it on top of the current time)
                systems_to_keep.append(failuretimes)  # keep track of these systems, they can still lose data if a replacement fails soon enough
            instances=instances+  sys_lost_data
            instances_wuer=instances_wuer+  sys_lost_data_wuer
    return   systems_to_keep,instances,instances_wuer

def MTTDL_calc(t_shards, p_shards, failure_rate,repair_time):
       d_shards=t_shards-p_shards
       cval=failure_rate*(failure_rate/repair_time)**p_shards /((d_shards)*factorial(t_shards)/(factorial(d_shards)*factorial(p_shards)));
      # print("MTTDL for "+str(t_shards) +" total shards with "+
       # str(p_shards)+" shards, i.e.,"+ str(d_shards) +"+"+str(p_shards)+ ",  "+str(failure_rate)+"% AFR,  recovery speed "+ "("+str(round(repair_time,1)) +" days). -->"
       # +str(cval)+" ~"+str(-math.log10(365/cval)) )
       return   cval


def MTTDL_calc_uer(t_shards, p_shards, failure_rate,repair_time,hval):
       MTTDL_th_cred=MTTDL_calc(t_shards, p_shards-1, failure_rate,repair_time)
       MTTDL_th_c=MTTDL_calc(t_shards, p_shards, failure_rate,repair_time)
       invMTTL= 1/MTTDL_th_c + hval/MTTDL_th_cred
       return 1/invMTTL;
def h_calc(d_shards,uer,d_cap):

       return 1-math.exp(-uer*d_shards*d_cap*8/1000)


def simulate(t_shards, p_shards, afr,uer, d_cap, r_speed,simulation_size_scale):
    d_shards =t_shards -p_shards
    repair_time=10**6*d_cap/r_speed/(60*60*24)
    failure_rate= -365.25/math.log(1-afr/100)     # this is in 1/days
    hval= h_calc(d_shards,uer,d_cap)
    MTTDL_th=MTTDL_calc(t_shards,p_shards, failure_rate,repair_time)
    MTTDL_th_uer=MTTDL_calc_uer(t_shards,p_shards, failure_rate,repair_time,hval)

    nines_th= -math.log10(365/MTTDL_th)
    nines_th_uer= -math.log10(365/MTTDL_th_uer)
    sim_size=simulation_size_scale*max(10000,min(50*10**math.ceil(nines_th) ,40000000) )
    if p_shards>1:
        ptext=" parities"
    else:
        ptext=" parity"
    print("Running "+str(sim_size)+" simulations: "+str(t_shards) +" total shards with "+
        str(p_shards)+ptext+", i.e.,"+ str(d_shards) +"+"+str(p_shards)+ ",  "+str(afr)+"% AFR, uer="+ str(uer)+" * 10^-{15}, "+str(d_cap)+"TB drive capacity, and "+
        str(r_speed)+"MB/s recovery speed "+ "("+str(round(repair_time,1)) +" days)." )
    print(  "prob of UER="+str(round(hval,2)   )    +     ". Theoretical predictions:  NoN= "+    str(round(nines_th,2))+" nines, NoN wUER= " +    str(round(nines_th_uer,2))+" nines. " )

    #initiate the full simulation data
    systems=[]
    for sysid in range(0,sim_size):
        systems.append(np.random.exponential(failure_rate, t_shards))

    totalinstances=0
    totalinstances_wuer=0
    while len(systems)>0:
        returned=sifter(systems,p_shards,repair_time,failure_rate,hval)
        systems=returned[0]
        totalinstances=totalinstances+returned[1]
        totalinstances_wuer=totalinstances_wuer+returned[2]

   # print("Simulated "+str(sim_size)+ " systems with " +str(totalinstances ))
    if totalinstances>0:
        print("Simulation Results:"+str(totalinstances )+"  data loss instances: NoN= "
        +str(round(-math.log10(totalinstances/sim_size),2))+" nines, NoN wUER= "+ str(round(-math.log10(totalinstances_wuer/sim_size),2))+ " nines" )
    else:
        print("No failures detected. Try to increase the simulation size!")

     #with one sigma="+str(round(1/(totalinstances**0.5),2))
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('number_of_total_shards', type=int),
    parser.add_argument('number_of_parity_shards', type=int),
    parser.add_argument('annual_failure_rate_in_pct', type=float),
    parser.add_argument('uer', type=float),   # in the units of 10^{-15}
    parser.add_argument('drive_capacity_in_TB', type=float),
    parser.add_argument('recovery_speed_in_MBps', type=float),
    parser.add_argument('simulation_size_scale', type=int),
    args = parser.parse_args()
    simulate(args.number_of_total_shards, args.number_of_parity_shards, args.annual_failure_rate_in_pct,args.uer, args.drive_capacity_in_TB, args.recovery_speed_in_MBps,args.simulation_size_scale)


if __name__ == '__main__':
    main()