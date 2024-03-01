import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.size"] = 20
plt.rcParams["axes.labelsize"] = 20
#label_color = "#595959"
label_color = "black"
plt.rcParams["text.color"] = label_color
plt.rcParams["axes.labelcolor"] = label_color
plt.rcParams['xtick.color'] = label_color
plt.rcParams['ytick.color'] = label_color

red_color = "#C00000"
blue_color = "#0070C0"
grey_color = "#A5A5A5"

gs_es_dic = pd.read_excel("gs_es.xlsx")

def plot_figures(output_dir, exclude=False, savefigs=True):
    exclusions = ["3ejh_A", "3m7p_A",
           "1cee_B", "2k42_A",
           "2kb8_A", "6vw2_A",
           "4cmq_B", "4zt0_C",
           "3t5o_A", "4a5w_B",
           "1x0g_A", "1x0g_D",
           "1wyy_B", "5wrg_C",
           "1g2c_F", "5tpn_A",
           "1svf_C", "4wsg_C",
           "1wp8_C", "5ejb_C",
           "4h2a_A", "3j9c_A",
           "3mko_A", "5ine_A",
           "2lv1_A", "6lni_A",
           "2kkw_A", "2n0a_D",
           "2mwf_A", "2nnt_A",
           "1iyt_A", "2nao_F"]

    if not exclude:
        exclusions = []

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    metrics = ["tm", "plddt", "ptm", "composite", "fsr_tm", "fsr_plddt"]
    metric_titles = {"tm": "Global TM", "plddt": "Global pLDDT", "ptm": "pTM", "composite": "composite", "fsr_tm": "FSR TM", "fsr_plddt": "FSR pLDDT"}

    for m in metrics:
        plt.figure(figsize=(9,6))

        x = [gs_es_dic[m+"_gs"][i] for i in range(len(gs_es_dic[m+"_gs"])) if gs_es_dic["ground_states"][i] != "equilibrium" and gs_es_dic["pdb_gs"][i] not in exclusions]
        y = [gs_es_dic[m+"_es"][i] for i in range(len(gs_es_dic[m+"_es"])) if gs_es_dic["ground_states"][i] != "equilibrium" and gs_es_dic["pdb_gs"][i] not in exclusions]
        if m == "fsr_tm":
            fsr_tm_saved_gs = x
            fsr_tm_saved_es = y
        if m == "fsr_plddt":
            composite_saved_gs = x
            composite_saved_es = y

        col = [blue_color if x[i0]>y[i0] else red_color for i0 in range(len(x))]
        x_equil = [gs_es_dic[m+"_gs"][i] for i in range(len(gs_es_dic[m+"_gs"])) if gs_es_dic["ground_states"][i] == "equilibrium" and gs_es_dic["pdb_gs"][i] not in exclusions]
        y_equil = [gs_es_dic[m+"_es"][i] for i in range(len(gs_es_dic[m+"_es"])) if gs_es_dic["ground_states"][i] == "equilibrium" and gs_es_dic["pdb_gs"][i] not in exclusions]

        if m in ["tm", "ptm", "fsr_tm"]:
            plt.plot([0,1],[0,1], "k")  
            plt.xlim([-0.05,1.05])
            plt.ylim([-0.05,1.05])
        else:
            plt.plot([0,100],[0,100], "k")  
            plt.xlim([-5,105])
            plt.ylim([-5,105])

        plt.scatter(x, y, s=60, c = col)
        plt.scatter(x_equil, y_equil, s=60, c = grey_color)
        plt.tick_params(top=True, right=True)
        #plt.title(metric_titles[m], fontweight="bold")
        plt.xlabel("Lower energy conformation")
        plt.ylabel("Higher energy conformation")
        if savefigs:
            plt.savefig(output_dir+"/"+m+".png", dpi=300) 
            plt.clf()

        if exclude:
            accuracy = sum([1 if x[i]>y[i] else 0 for i in range(len(x))])/len(x)
            print("{:} accuracy: {:}".format(m, accuracy))
    '''plt.scatter(fsr_tm_saved_gs,composite_saved_gs)
    plt.xlabel("fsr_tm")
    plt.ylabel("composite")
    plt.savefig(output_dir+"/fsr_plddt_vs_tm_gs.png")
    plt.clf()
    plt.scatter(fsr_tm_saved_es,composite_saved_es)
    plt.xlabel("fsr_tm")
    plt.ylabel("composite")
    plt.savefig(output_dir+"/fsr_plddt_vs_tm_es.png")'''

def main():
    plot_figures("figures")
    plot_figures("figures_exclusions", exclude=True)


main()