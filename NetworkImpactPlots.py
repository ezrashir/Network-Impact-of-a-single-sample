from steadystates import *
import scipy.stats as scs
from scipy.spatial.distance import jensenshannon
import pandas as pd
from sklearn import metrics


def create_net_intersected(mat:np.array, sig_p_val:float=1e-3):
    mat_size = np.size(mat[0, :])
    network = np.zeros((mat_size, mat_size))
    for i in range(mat_size):
        j_val = np.arange(i, mat_size, 1)
        for j in j_val:
            intersect_indexes = np.intersect1d(np.nonzero(mat[:, i]), np.nonzero(mat[:, j]))
            vec1, vec2 = [], []
            for k in range(np.size(intersect_indexes)):
                vec1.append(mat[intersect_indexes[k], i])
                vec2.append(mat[intersect_indexes[k], j])

            x, p = scs.pearsonr(vec1, vec2)
            if p < sig_p_val:
                network[i, j] = x
    network = network - (np.identity(mat_size))
    return network
def Network_Impact_parameters_Calculate(reference_cohort, test_set):
    # test_set should contain all samples composing test network (i.e. reference cohort + test sample)

    def binaryTrans(vec):
        for i in range(np.size(vec)):
            if vec[i] != 0:
                vec[i] = 1
        return vec

    inter_val, dropped_links, unchanged_links = [], 0, 0
    reference_flat = np.matrix.flatten(create_net_intersected(reference_cohort))
    test_flat = np.matrix.flatten(create_net_intersected(test_set))

    intersect_indexes = np.intersect1d(np.nonzero(reference_flat), np.nonzero(test_flat))
    for j in intersect_indexes:
        inter_val0 = reference_flat[j] - test_flat[j]
        inter_val.append(inter_val0)

        if np.abs(reference_flat[j]) > np.abs(test_flat[j]):
            dropped_links += 1
        if np.abs(reference_flat[j]) == np.abs(test_flat[j]):
            unchanged_links += 1

    DW = np.mean(inter_val)
    DS = 1 - metrics.jaccard_score(binaryTrans(reference_flat), binaryTrans(test_flat))
    changed_common_links = np.size(intersect_indexes) - unchanged_links
    theta = dropped_links / changed_common_links
    return DS, DW, theta


def semi_supervised_plots(file_path:str):
    reference_cohort = pd.read_excel(file_path, sheet_name='Reference cohort', header=None).to_numpy()
    GLV_SS_cohort = pd.read_excel(file_path, sheet_name='GLV Steady States cohort', header=None).to_numpy()
    shuffled_cohort = pd.read_excel(file_path, sheet_name='Shuffled cohort', header=None).to_numpy()

    num_of_test_samples = np.size(GLV_SS_cohort, 0)

    def calculate_theta_semi_supervised():
        theta_SS, theta_sh = np.zeros(num_of_test_samples), np.zeros(num_of_test_samples)
        for i in range(num_of_test_samples):
            _, _, theta_SS[i] = Network_Impact_parameters_Calculate(reference_cohort=reference_cohort,
                                                                    test_set=np.vstack(
                                                                        (reference_cohort, GLV_SS_cohort[i, :])))
            _, _, theta_sh[i] = Network_Impact_parameters_Calculate(reference_cohort=reference_cohort,
                                                                    test_set=np.vstack(
                                                                        (reference_cohort, shuffled_cohort[i, :])))
        return theta_SS, theta_sh

    def calculate_rJSD():
        rJSD_SS, rJSD_sh = np.zeros(num_of_test_samples), np.zeros(num_of_test_samples)
        for i in range(num_of_test_samples):
            dist_GLV, dist_sh = [], []
            for j in range(np.size(reference_cohort, 0)):
                dist_GLV.append(jensenshannon(GLV_SS_cohort[i, :], reference_cohort[j, :]))
                dist_sh.append(jensenshannon(shuffled_cohort[i, :], reference_cohort[j, :]))
            rJSD_SS[i] = np.mean(dist_GLV)
            rJSD_sh[i] = np.mean(dist_sh)
        return rJSD_SS, rJSD_sh

    def ROC_compare(dist0_positive:np.array, dist0_negative:np.array, dist1_positive:np.array, dist1_negative:np.aray, threshold_span:float):
        def ROC_data(dist_positive:np.array, dist_negative:np.array, threshold_span:float, bins_number:int=1000):
            fpr, tpr = [], []
            thresholds = np.arange(threshold_span[0], threshold_span[1],
                                   (threshold_span[1] - threshold_span[0]) / bins_number)
            for i in thresholds:
                fpr.append(np.size(np.where(dist_positive > i)) / np.size(dist_positive))
                tpr.append(np.size(np.where(dist_negative > i)) / np.size(dist_positive))
            auc = metrics.auc(fpr, tpr)
            return fpr, tpr, auc

        # dist0 (positive and negative, or GLV_SS and shuffled profiles) are the theta distributions and dist1 are the two rJSD distributions
        fpr0, tpr0, auc0 = ROC_data(dist0_positive, dist0_negative, threshold_span=threshold_span)
        fpr1, tpr1, auc1 = ROC_data(dist1_positive, dist1_negative, threshold_span=threshold_span)
        return fpr0, tpr0, auc0, fpr1, tpr1, auc1

    theta_SS, theta_sh = calculate_theta_semi_supervised()
    rJSD_SS, rJSD_sh = calculate_rJSD()
    fpr_SS, tpr_SS, AUC_SS, fpr_sh, tpr_sh, AUC_sh = ROC_compare(theta_SS, theta_sh, rJSD_SS, rJSD_sh,
                                                                 threshold_span=[0, 1])

    fig, ax = plt.subplots(1, 3)
    plt.subplots_adjust(left=0.05, bottom=None, right=0.96, top=None, wspace=None, hspace=None)
    font_size = 10

    ax[0].hist(rJSD_SS, histtype='stepfilled', alpha=0.5, label='GLV', color='darkorange')
    ax[0].hist(rJSD_sh, histtype='stepfilled', alpha=0.5, label='shuffled', color='olivedrab')
    ax[0].legend()
    ax[0].set_title('rJSD values')
    ax[1].hist(theta_SS, histtype='stepfilled', alpha=0.5, label='GLV', color='darkorange')
    ax[1].hist(theta_sh, histtype='stepfilled', alpha=0.5, label='shuffled', color='olivedrab')
    ax[1].legend()
    ax[1].set_title(r'$\theta$ values')
    ax[2].plot(fpr_SS, tpr_SS, color='lightseagreen', label=r'$\theta$', linewidth=3)
    ax[2].plot(fpr_sh, tpr_sh, color='magenta', label='rJSD', linewidth=3)
    ax[2].legend()
    ax[2].set_title('ROC curves')
    ax[2].set_xlabel('False positive rate', fontsize=font_size)
    ax[2].set_ylabel('True positive rate', fontsize=font_size)
    ax[2].text(0.65, 0.5, r'AUC = ' + str(np.round(AUC_SS, decimals=2)), color='darkcyan', fontsize=font_size)
    ax[2].text(0.65, 0.4, r'AUC = ' + str(np.round(AUC_sh, decimals=2)), color='magenta', fontsize=font_size)
    plt.show()
def supervised_plots(file_path):
    reference_cohort_A = pd.read_excel(file_path, sheet_name='Reference cohort A', header=None).to_numpy()
    test_cohort_A = pd.read_excel(file_path, sheet_name='Test cohort A', header=None).to_numpy()
    reference_cohort_B = pd.read_excel(file_path, sheet_name='Reference cohort B', header=None).to_numpy()
    test_cohort_B = pd.read_excel(file_path, sheet_name='Test cohort B', header=None).to_numpy()

    m = np.size(reference_cohort_A, 0)
    test_size = np.size(test_cohort_A, 0)

    def calculate_network_impact_supervised(reference_cohort:np.array, test_samples:np.array):
        DS, DW, theta = [], [], []
        for i in range(np.size(test_samples, 0)):
            test_set = np.vstack([reference_cohort, test_samples[i, :]])
            DS0, DW0, theta0 = Network_Impact_parameters_Calculate(reference_cohort, test_set)
            DS.append(DS0)
            DW.append(DW0)
            theta.append(theta0)
        return DS, DW, theta

    def calculate_rJSD_supervised(reference_cohort:np.array, test_samples:np.array):
        rJSD = np.zeros(test_size)
        k = 0
        for j in range(test_size):
            rJSD0 = []
            for i in range(m):
                rJSD0.append(jensenshannon(reference_cohort[i, :], test_samples[j, :]))
            rJSD[k] = np.mean(rJSD0)
            k += 1
        return rJSD

    DS_A2A, DW_A2A, theta_A2A = calculate_network_impact_supervised(reference_cohort_A, test_cohort_A)
    DS_A2B, DW_A2B, theta_A2B = calculate_network_impact_supervised(reference_cohort_A, test_cohort_B)
    DS_B2A, DW_B2A, theta_B2A = calculate_network_impact_supervised(reference_cohort_B, test_cohort_A)
    DS_B2B, DW_B2B, theta_B2B = calculate_network_impact_supervised(reference_cohort_B, test_cohort_B)

    rJSD_A2A = calculate_rJSD_supervised(reference_cohort_A, test_cohort_A)
    rJSD_A2B = calculate_rJSD_supervised(reference_cohort_A, test_cohort_B)
    rJSD_B2A = calculate_rJSD_supervised(reference_cohort_B, test_cohort_A)
    rJSD_B2B = calculate_rJSD_supervised(reference_cohort_B, test_cohort_B)

    fig, ax = plt.subplots(2, 2)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.35, hspace=None)

    ax[0, 0].scatter(rJSD_A2A, rJSD_B2A, label='GLV model A', color='darkorange')
    ax[0, 0].scatter(rJSD_A2B, rJSD_B2B, label='GLV model B', color='lightslategrey')
    ax[0, 0].set_xlabel('rJSD w.r.t model A')
    ax[0, 0].set_ylabel('rJSD w.r.t model B')
    ax[0, 0].set_title('rJSD')
    ax[0, 0].legend()

    ax[0, 1].scatter(theta_A2A, theta_B2A, label='d for tests from dynamics 1\A', color='darkorange')
    ax[0, 1].scatter(theta_A2B, theta_B2B, label='d for tests from dynamics 2\B', color='lightslategrey')
    ax[0, 1].set_xlabel(r'$\theta$ w.r.t model A')
    ax[0, 1].set_ylabel(r'$\theta$ w.r.t model B')
    ax[0, 1].set_title(r'$\theta$')

    ax[1, 0].scatter(DS_A2A, DS_B2A, label='a for tests from dynamics 1\A', color='darkorange')
    ax[1, 0].scatter(DS_A2B, DS_B2B, label='a for tests from dynamics 2\B', color='lightslategrey')
    ax[1, 0].set_xlabel(r'$\Delta$S w.r.t model A')
    ax[1, 0].set_ylabel(r'$\Delta$S w.r.t model B')
    ax[1, 0].set_title(r'$\Delta$S')

    ax[1, 1].scatter(DW_A2A, DW_B2A, label='d for tests from dynamics 1\A', color='darkorange')
    ax[1, 1].scatter(DW_A2B, DW_B2B, label='d for tests from dynamics 2\B', color='lightslategrey')
    ax[1, 1].ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
    ax[1, 1].set_xlabel(r'$\Delta$W w.r.t model A')
    ax[1, 1].set_ylabel(r'$\Delta$W w.r.t model B')
    ax[1, 1].set_title(r'$\Delta$W')

    plt.show()


# To plot the semi supervised Fig 3.(a-c) insert the file path on your computer to the following function:
# The plot function takes about 6-8 minutes from file to figure for the cohorts in the provided .xlsx file
# in which reference cohort size: m=100 and there are 100 test samples of each kind (GLV steady state and 
# shuffled profiles).
semi_supervised_plots('Insert The Path To The File From Your Computer')


# To plot the supervised Fig 5.(a-d) insert the file path on your computer to the following function:
# The plot function takes about 10-15 minutes from file to figure for the cohorts in the provided .xlsx file
# in which reference cohort size (in both models) is m=50 and there are 100 test samples of each kind (GLV model A and GLV model B).
supervised_plots('Insert The Path To The File From Your Computer')


