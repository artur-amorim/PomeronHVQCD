import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

def plotF2(F2_pred, F2_exp_data, output_file):
    """
        Plots the F2 predictions vs F2 data and saves the result
        in a pdf file called output_file.pdf
        The plot is splitted in two subplots with Q2 < 8.5 and Q2 >= 8.5.
        F2_pred : pandas data frame that contains the predictions of F2.
        F2_exp_data : pandas data frame that contains the experimental data of F2
        output_file : string that stores the name of output file 
    """
    # Split the data frame
    F2_pred_Q2max_85 = F2_pred[F2_pred.Q2 < 8.5]
    F2_pred_Q2min_85 = F2_pred[F2_pred.Q2 >= 8.5]
    # Plot the subplot with Q2 < 8.5
    fig, ax = plt.subplots(2,1,sharex=False, sharey=False, figsize = (7,10))
    x_min = 10**(-6)
    x_max = 0.01
    ax[0].set_xlim([x_min, x_max])
    ax[0].set_xlabel("x")
    ax[0].set_xscale("log")
    colours = ["r","b","g"]
    Q2s = F2_pred_Q2max_85.Q2.unique()
    n = len(Q2s)
    for i in range(n) :
        # plot the predicted points
        x_values = np.array(F2_pred_Q2max_85[F2_pred_Q2max_85.Q2 == Q2s[i]].x)
        struct_func_values = np.array(F2_pred_Q2max_85[F2_pred_Q2max_85.Q2 == Q2s[i]].F2)
        ax[0].plot(x_values, struct_func_values, colours[i%3], linewidth = 0.5)
        # plot the experimental points
        x_values = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].x)
        struct_func_values = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].F2)
        struct_func_errors = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].F2err)
        pt_fmt = colours[i%3] + "o"
        ax[0].errorbar(x_values, struct_func_values, struct_func_errors, fmt=pt_fmt, markersize=1)
        x_annotate = 0.5 * x_values.min()
        y_annotate = 1.1 * struct_func_values.max()
        text = ax[0].annotate(str(Q2s[i]), xy = (x_annotate, y_annotate),
                                xytext=(x_annotate, y_annotate), color = colours[i%3])
        text.set_fontsize(6)
    ax[0].set_ylim([0.1, 1.4])
    x_min = 1.5*10**(-4)
    ax[1].set_xlim([x_min, x_max])
    ax[1].set_xlabel("x")
    ax[1].set_xscale("log")
    Q2s = F2_pred_Q2min_85.Q2.unique()
    n = len(Q2s)
    for i in range(n) :
        # plot the predicted points
        x_values = np.array(F2_pred_Q2min_85[(F2_pred_Q2min_85.Q2 == Q2s[i]) & (F2_pred_Q2min_85.x >= x_min)].x)
        struct_func_values = np.array(F2_pred_Q2min_85[(F2_pred_Q2min_85.Q2 == Q2s[i]) & (F2_pred_Q2min_85.x >= x_min)].F2)
        ax[1].plot(x_values, struct_func_values, colours[i%3], linewidth = 0.5)
        # plot the experimental points
        x_values = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].x)
        struct_func_values = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].F2)
        struct_func_errors = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].F2err)
        pt_fmt = colours[i%3] + "o"
        ax[1].errorbar(x_values, struct_func_values, struct_func_errors, fmt=pt_fmt, markersize=1)
        x_annotate = 2.0 * 10**(-3)
        y_annotate = np.min(F2_pred_Q2min_85[(F2_pred_Q2min_85.Q2 == Q2s[i]) & (F2_pred_Q2min_85.x - x_annotate < 1.5*10**(-4))].F2)
        text = ax[1].annotate(str(Q2s[i]), xy = (x_annotate, y_annotate),
                                                    xytext=(x_annotate, y_annotate), color = colours[i%3])
        text.set_fontsize(6)
    ax[1].set_ylim(0.6, 1.75)
    # Some cosmetic surgery
    fig.tight_layout()
    fig.suptitle(r'$F_2(x,Q^2)$ vs x')
    fig.subplots_adjust(top=0.95)
    # Save the plot
    fig.savefig(output_file+".pdf", dpi = 100)


print("Program usage: python3 plotF2.py pred_data exp_data output_file_name")
arguments = sys.argv
# Load experimental and predicted data
F2_pred = pd.read_csv(arguments[1], sep = '\t')
F2_exp = pd.read_csv(arguments[2], sep = '\t')
plotF2(F2_pred, F2_exp, arguments[3])