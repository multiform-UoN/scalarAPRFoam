#!/usr/bin/python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import signal
from scipy.integrate import solve_ivp
from lmfit import Parameters, Minimizer, report_fit

# -----------------------------
# Read Excel data
# -----------------------------
def read_excel_data(file_path):
    xls = pd.ExcelFile(file_path, engine='openpyxl')
    sheets_data = {}
    species_names = {}
    sheet_params = {}

    for sheet_name in xls.sheet_names:
        data = pd.read_excel(xls, sheet_name=sheet_name)
        time = data.iloc[:, 0].values
        concentrations = data.iloc[:, 1:5].values
        column_names = data.columns[1:]
        species_names[sheet_name] = [col for col in column_names if not (str(col).startswith('Unnamed') or str(col).startswith('DATA'))]

        # ---- group by time (mean) ----
        df = pd.DataFrame(concentrations, columns=species_names[sheet_name])
        df['time'] = time
        grouped = df.groupby('time')
        time_unique = grouped['time'].first().values
        mean_conc = grouped.mean().drop(columns=['time'], errors='ignore').values
        std_conc = grouped.std().drop(columns=['time'], errors='ignore').fillna(0).values

        sheets_data[sheet_name] = (time_unique, mean_conc, std_conc)

        T = data.iloc[0, 7]
        Vliq = data.iloc[1, 7] * 1e-6
        Vgas = data.iloc[2, 7] * 1e-6
        kgcat = data.iloc[3, 7] * 1e-3
        sheet_params[sheet_name] = (T, Vliq, Vgas, kgcat)

    return sheets_data, species_names, sheet_params

# -----------------------------
# Reaction model
# -----------------------------
def reaction_modelLOG(t, y, kineticParameters, T, Vliq, Vgas, kgcat):
    G, P, H, C = y
    k_APR, k_pg, Kg, Kh = kineticParameters
    R = 8.3144

    henry = 7.80e-6 * np.exp(530 * (1 / T - 1 / 298))
    henry = henry * R * T

    k_APR = 10**k_APR
    k_pg = 10**k_pg
    Kg = 375e-3
    Kh = 10**Kh

    den = 1 + Kg * G + np.sqrt(Kh * H * henry)  if H > 0 else 1 + Kg * G

    r1 = k_APR * G / den**2 if G > 0 else 0
    r2 = k_pg * G * ((H * henry)**0.5) / den**2 if (G > 0 or H > 0) else 0

    dG_dt = (kgcat/Vliq) *  (-r1 - r2)
    dP_dt = (kgcat/Vliq) *  r2
    dH_dt = (kgcat/Vliq) *  (7 * r1 - r2) * Vliq / Vgas
    dC_dt = (kgcat/Vliq) *  (3 * r1) * Vliq / Vgas

    return [dG_dt, dP_dt, dH_dt, dC_dt]

def reaction_modelLIN(t, y, kineticParameters, T, Vliq, Vgas, kgcat):
    G, P, H, C = y
    k_APR, k_pg, Kg, Kh = kineticParameters
    R = 8.3144

    henry = 7.80e-6 * np.exp(530 * (1 / T - 1 / 298))
    henry = henry * R * T

    den = 1 + Kg * G + np.sqrt(Kh * H * henry)  if H > 0 else 1 + Kg * G

    r1 = k_APR * G / den**2 if G > 0 else 0
    r2 = k_pg * G * ((H * henry)**0.5) / den**2 if (G > 0 or H > 0) else 0

    dG_dt = (kgcat/Vliq) *  (-r1 - r2)
    dP_dt = (kgcat/Vliq) *  r2
    dH_dt = (kgcat/Vliq) *  (7 * r1 - r2) * Vliq / Vgas
    dC_dt = (kgcat/Vliq) *  (3 * r1) * Vliq / Vgas

    return [dG_dt, dP_dt, dH_dt, dC_dt]
# -----------------------------
# ODE solver with timeout
# -----------------------------
def solve_ode(model, kineticParameters, initial_conditions, time_points, T, Vliq, Vgas, kgcat):
    class TimeoutException(Exception):
        pass

    def timeout_handler(signum, frame):
        raise TimeoutException("solve_ivp exceeded the time limit")

    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(60)

    try:
        solution = solve_ivp(
            model,
            [time_points[0], time_points[-1]],
            initial_conditions,
            t_eval=time_points,
            args=(kineticParameters, T, Vliq, Vgas, kgcat),
            method='BDF',
            atol=1e-8,
            rtol=1e-8,
            max_step=60,
            # min_step=1e-12
        )
    except (TimeoutException, ValueError) as e:
        print(f"solve_ivp failed: {e}")
        return np.full((len(time_points), len(initial_conditions)), np.nan)
    finally:
        signal.alarm(0)

    if solution.status < 0 or not solution.success:
        print(f"Integration failed: {solution.message}")
        return np.full((len(time_points), len(initial_conditions)), np.nan)

    return solution.y.T

# -----------------------------
# Residual function for lmfit
# -----------------------------
def residualsASS(params, sheets_data, initial_conditions, sheet_params, model):
    kineticParameters = [
        params['k_APR'].value,
        params['k_pg'].value,
        params['Kg'].value,
        params['Kh'].value
    ]

    resids = []
    for sheet_name, (time_exp, conc_exp, std_exp) in sheets_data.items():
        T, Vliq, Vgas, kgcat = sheet_params[sheet_name]
        conc_model = solve_ode(model, kineticParameters, initial_conditions[sheet_name], time_exp, T, Vliq, Vgas, kgcat)

        if np.isnan(conc_model).any():
            # Return large residuals with SAME SHAPE as expected data
            large_resids = np.ones_like(conc_exp) * 1e12
            resids.append(large_resids.flatten())
        else:
            # Use relative residuals with small regularization to avoid division by zero
            relative_resids = np.abs((conc_exp - conc_model))# / (conc_exp + 1e-8))
            resids.append(relative_resids.flatten())

    all_resids = np.concatenate(resids)
    print(f"Current params: {params['k_APR'].value:.2e} {params['k_pg'].value:.2e} {params['Kg'].value:.2e} {params['Kh'].value:.2e} --> SSE: {np.sum(all_resids**2):.3e}")
    return all_resids

def residualsREL(params, sheets_data, initial_conditions, sheet_params, model):
    kineticParameters = [
        params['k_APR'].value,
        params['k_pg'].value,
        params['Kg'].value,
        params['Kh'].value
    ]

    resids = []
    for sheet_name, (time_exp, conc_exp, std_exp) in sheets_data.items():
        T, Vliq, Vgas, kgcat = sheet_params[sheet_name]
        conc_model = solve_ode(model, kineticParameters, initial_conditions[sheet_name], time_exp, T, Vliq, Vgas, kgcat)

        if np.isnan(conc_model).any():
            # Return large residuals with SAME SHAPE as expected data
            large_resids = np.ones_like(conc_exp) * 1e12
            resids.append(large_resids.flatten())
        else:
            # Use relative residuals with small regularization to avoid division by zero
            relative_resids = ((conc_exp - conc_model) / (conc_exp + 1e-8))**2
            resids.append(relative_resids.flatten())

    all_resids = np.sum(np.concatenate(resids))
    print(f"Current params: {params['k_APR'].value:.2e} {params['k_pg'].value:.2e} {params['Kg'].value:.2e} {params['Kh'].value:.2e} --> RelErr: {all_resids:.3e}")
    return all_resids

# -----------------------------
# Plot results with means + std
# -----------------------------
def plot_results(sheets_data, species_names, result, initial_conditions, sheet_params):
    kineticParameters = [
        result.params['k_APR'].value,
        result.params['k_pg'].value,
        result.params['Kg'].value,
        result.params['Kh'].value
    ]

    all_exp, all_exp_err, all_model, all_colors, all_shapes = [], [], [], [], []
    markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x', '|', '_']
    colors = plt.cm.tab10.colors
    legend_elements = []

    for sheet_idx, (sheet_name, (time_exp, conc_exp, std_exp)) in enumerate(sheets_data.items()):
        T, Vliq, Vgas, kgcat = sheet_params[sheet_name]

        # model at dense grid
        conc_model = solve_ode(reaction_modelLIN, kineticParameters,initial_conditions[sheet_name],np.linspace(time_exp[0], time_exp[-1], 100), T, Vliq, Vgas, kgcat)
        species_labels = species_names[sheet_name]

        # Create figure with dual y-axes
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax2 = ax1.twinx()  # Create secondary y-axis
        
        # Assume glycerol is the first species (index 0)
        glycerol_idx = 0
        
        for i in range(conc_exp.shape[1]):
            if i == glycerol_idx:  # Plot glycerol on secondary axis
                # Experimental data for glycerol
                ax2.errorbar(
                    time_exp/3600,
                    conc_exp[:, i],
                    yerr=std_exp[:, i],
                    fmt='x',
                    color=colors[i % len(colors)],
                    markersize=5,
                    markerfacecolor='none',
                    label=f"Exp {species_labels[i]}"
                )

                # Model curve for glycerol
                ax2.plot(
                    np.linspace(time_exp[0], time_exp[-1], 100)/3600, 
                    conc_model[:, i], 
                    linestyle='-', 
                    color=colors[i % len(colors)],
                    label=f"Model {species_labels[i]}"
                )
                
                # Set labels and styling for secondary axis
                ax2.set_ylabel("Glycerol Concentration, mol/m³", color=colors[glycerol_idx])
                ax2.tick_params(axis='y', labelcolor=colors[glycerol_idx])
                
            else:  # Plot other species on primary axis
                # Experimental data for other species
                ax1.errorbar(
                    time_exp/3600,
                    conc_exp[:, i],
                    yerr=std_exp[:, i],
                    fmt='x',
                    color=colors[i % len(colors)],
                    markersize=5,
                    markerfacecolor='none',
                    label=f"Exp {species_labels[i]}"
                )

                # Model curve for other species
                ax1.plot(
                    np.linspace(time_exp[0], time_exp[-1], 100)/3600, 
                    conc_model[:, i], 
                    linestyle='-', 
                    color=colors[i % len(colors)],
                    label=f"Model {species_labels[i]}"
                )

            # parity plot data (with error bars) - unchanged
            model_interp = np.interp(time_exp, np.linspace(time_exp[0], time_exp[-1], 100), conc_model[:, i])

            all_exp.extend(conc_exp[:, i])
            all_exp_err.extend(std_exp[:, i])
            all_model.extend(model_interp)
            all_colors.extend([colors[sheet_idx % len(colors)]] * len(conc_exp[:, i]))
            all_shapes.extend([markers[i % len(markers)]] * len(conc_exp[:, i]))

            legend_elements.append(plt.Line2D([0], [0], color=colors[sheet_idx % len(colors)], marker=markers[i % len(markers)], linestyle='None', label=f"{sheet_name} - {species_labels[i]}"))

        # Set labels for primary axis
        ax1.set_xlabel("Time, h")
        ax1.set_ylabel("Concentration, mol/m³")
        ax1.set_title(f"Experimental vs Model Data - {sheet_name}")
        
        # Combine legends from both axes
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')
        
        plt.tight_layout()
        plt.savefig("plot_" + sheet_name + ".png", dpi=300, transparent=True, bbox_inches='tight')
        plt.close(fig)  # Close figure to free memory

    # parity plot - unchanged
    plt.figure("Parity Plot")
    for exp, exp_err, model, color, shape in zip(all_exp, all_exp_err, all_model, all_colors, all_shapes):
        plt.errorbar(exp, model, xerr=exp_err, fmt=shape, color=color, alpha=0.7)
    plt.plot([min(all_exp), max(all_exp)], [min(all_exp), max(all_exp)], 'r--', label="Parity Line")
    plt.xlabel("Experimental Concentration, mol/m^3")
    plt.ylabel("Model Concentration, mol/m^3")
    plt.title("Parity Plot - Experimental vs Model")
    plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', title="Legend")
    plt.savefig("parity_plot.png", dpi=300, transparent=True, bbox_inches='tight')

# -----------------------------
# Save results
# -----------------------------
def write_results_lmfit(result):
    with open("output_file.txt", 'w') as f:
        for name, par in result.params.items():
            f.write(f"{name}: {par.value:.6e} ± {par.stderr if par.stderr else np.nan:.6e}\n")
        f.write(f"Kg: 3.75e-1 ± 0.25\n")
        f.write(f"Chi-square: {result.chisqr}\n")
        f.write(f"Reduced Chi-square: {result.redchi}\n")

# -----------------------------
# Main
# -----------------------------
if __name__ == "__main__":
    file_path = "rawDataComplite.xlsx"
    sheets_data, species_names, sheet_params = read_excel_data(file_path)
    if not sheets_data:
        raise SystemExit("No sheets found in the Excel file")

    initial_conditions = {sheet: data[1][0] for sheet, data in sheets_data.items()}

    initial_population = [[3.2, 3.4, 0.375, 8.7]]  # k_APR, k_pg, Kg, Kh, Kwgs

    params = Parameters()
    params.add('k_APR', value=initial_population[0][0], min=-3, max=6)
    params.add('k_pg', value=initial_population[0][1], min=-3, max=6)
    params.add('Kg', value=initial_population[0][2], min=0.3749999, max=0.3750001)
    params.add('Kh', value=initial_population[0][3], min=0, max=15)


    mini = Minimizer(residualsASS, params, fcn_args=(sheets_data, initial_conditions, sheet_params, reaction_modelLOG))

    global_result = mini.minimize(
        method='differential_evolution',
        params=params,
        popsize=15,
        mutation=(0.75,1.5),
        recombination=0.7,
        max_nfev=100000,
        strategy='best1bin',
        tol=1e-6,
        polish=False,
        disp=True,
        init='latinhypercube'
        )

    print("Global optimization finished. Best-fit parameters:")
    report_fit(global_result)

    params['k_APR'].set(value=10**global_result.params['k_APR'].value, min=1e-3, max=1e6)
    params['k_pg'].set(value=10**global_result.params['k_pg'].value, min=1e-3, max=1e6)
    params['Kg'].set(value=initial_population[0][2], min=0.275, max=0.475)
    params['Kh'].set(value=10**global_result.params['Kh'].value, min=1e-1, max=1e15)
    print("Refined parameters for local optimization:")
    report_fit(params)

    miniCI = Minimizer(residualsASS, params, fcn_args=(sheets_data, initial_conditions, sheet_params, reaction_modelLIN))
    # Use Levenberg-Marquardt (leastsq) method with enhanced stability
    result = miniCI.minimize(
        method = 'leastsq',
        params = params,  # Start from the best result of the global search
        # max_nfev = 10000,      # Maximum function evaluations (more iterations)
        # factor = 0.75,            # Initial step bound factor (smaller = more conservative)
        diag = None,            # Use automatic scaling
        ftol = 1e-12,            # Tolerance for termination by the change of the cost function
        xtol = 1e-12,            # Tolerance for termination by the change of the independent variables
        gtol = 1e-12             # Tolerance for termination by the norm of the gradient
    )

    print("Fit finished. Best-fit parameters:")
    report_fit(result)

    write_results_lmfit(result)
    plot_results(sheets_data, species_names, result, initial_conditions, sheet_params)

    print("Done. Plots saved and output_file.txt created.")
