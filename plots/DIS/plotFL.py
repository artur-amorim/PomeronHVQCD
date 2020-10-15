import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

def plotFL(FL_pred, FL_exp, output_file):
    """
        Plots the FL predictions vs FL data and saves the result
        in a pdf file called output_file.pdf
        The plot is splitted in two subplots with Q2 <= 10 and Q2 > 10.
        FL_pred : pandas data frame that contains the predictions of FL.
        FL_exp : pandas data frame that contains the experimental data of FL
        output_file : string that stores the name of output file 
    """
    # Split the FL_pred data frame
    Q2s_vals = FL_pred.Q2.unique()
    #FL_pred_Q2max_10 = FL_pred[FL_pred.Q2 <= 10]
    #FL_pred_Q2min_10 = FL_pred[FL_pred.Q2 > 10]
    # Find unique values of Q^2 in FL data with Q2 <= 10
    #Q2s_max_10 = FL_pred_Q2max_10.Q2.unique()
    #n_Q2max_10 = len(Q2s_max_10)
    # Find unique values of Q^2 in FL data with Q2 > 10
    #Q2s_min_10 = FL_pred_Q2min_10.Q2.unique()
    #n_Q2min_10 = len(Q2s_min_10)
    # number of rows
    n_rows = len(Q2s_vals)
    # Start the plot
    fig, ax = plt.subplots(n_rows,1,sharex=False, sharey=False, figsize = (7,12))
    fig.suptitle(r'$F_L(x,Q^2)$ vs x')
    #j_Q2max_10, j_Q2min_10 = 0, 0
    for i in range(n_rows) :
        figp = plt.subplot(n_rows,1,i+1)
        x_min, x_max = 0, 0
        if (Q2s_vals[i] < 12):
            x_min = 2*10**(-5)
            x_max = 8*10**(-4)
        else:
            x_min = 10**(-4)
            x_max = 3*10**(-3)
        figp.set_xlabel("x")
        x_values = np.array(FL_pred[FL_pred.Q2 == Q2s_vals[i]].x)
        struct_func_values = np.array(FL_pred[FL_pred.Q2 == Q2s_vals[i]].FL)
        plt.plot(x_values, struct_func_values, "b", linewidth = 0.5)
        plt.plot(x_values,0*x_values,"k--")
        x_annotate = 4*10**(-4)
        y_annotate = 0.5*(struct_func_values.min()+struct_func_values.max())
        text = plt.annotate(r'$Q^2 = ' + str(Q2s_vals[i]) + '$', xy = (x_annotate, y_annotate),
                                                    xytext=(1.1 * x_annotate, 0.7 * y_annotate))
        text.set_fontsize(10)
        # plot the experimental points
        x_values = np.array(FL_exp[FL_exp.Q2 == Q2s_vals[i]].x)
        struct_func_values = np.array(FL_exp[FL_exp.Q2 == Q2s_vals[i]].FL)
        struct_func_errors = np.array(FL_exp[FL_exp.Q2 == Q2s_vals[i]].FLerr)
        plt.errorbar(x_values, struct_func_values, struct_func_errors, fmt="bo", markersize=1)
        plt.xscale("log")
        plt.xlim(x_min, x_max)
    # Cosmetics
    l  = 0.1  # the left side of the subplots of the figure
    r = 0.975    # the right side of the subplots of the figure
    b = 0.05   # the bottom of the subplots of the figure
    t = 0.95      # the top of the subplots of the figure
    ws = 0.1   # the amount of width reserved for blank space between subplots
    hs = 0.25  # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left = l, right = r, top = t, bottom = b, wspace = ws, hspace = hs)
    plt.savefig(output_file+".pdf", dpi = 100)
    plt.close()



print("Program usage: python3 plotFL.py pred_data exp_data output_file_name")
arguments = sys.argv
# Load experimental and predicted data
FL_pred = pd.read_csv(arguments[1], sep = '\t')
FL_exp = pd.read_csv(arguments[2], sep = '\t')
plotFL(FL_pred, FL_exp, arguments[3])
