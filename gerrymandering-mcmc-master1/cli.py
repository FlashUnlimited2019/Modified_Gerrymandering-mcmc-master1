import argparse
from src.GerrymanderingMCMC import GerrymanderingMCMC
from mpi4py import MPI
import pandas as pd
import json
#default_file = "./src/data/prep4"
default_file = "./src/data/prep2_va"


def main():
    parser = argparse.ArgumentParser(
        description='Use MCMC Simulation to generate districting plans and plot relevant key statistics to illustrate the possibility that a source plan was gerrymandered')
    parser.add_argument("-g", "--graph_file", default=default_file,
                        help="A path to a potential districting plan specified in this projects proprietary json schema; defaults to ./src/data/iowa.json")
    parser.add_argument("-c", "--cooling_period", type=int, default=50,
                        help="The number of plans you'd like to generate _before_ counting them towards your ensemble; defaults to 50")
    parser.add_argument("-r", "--rounds", type=int, default=200,
                        help="The number of plans you'd like to generate and include in your ensemble; defaults to 200")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Include this flag if you'd like real-time output to the console")
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    graph_file = args.graph_file
    cooling_period = args.cooling_period
    rounds = args.rounds
    verbose = args.verbose
    cr = cooling_period+rounds
    # Build the gerrymandering MCMC using the variables you've been provided
    output = []
    for i in range(cr):
        if i == 0:
            mcmc = GerrymanderingMCMC(
                graph_file, cooling_period=cooling_period, verbose=verbose)
        else:
            mcmc = GerrymanderingMCMC(
                graph_file+str(rank), cooling_period=cooling_period, verbose=verbose)
        # mcmc = GerrymanderingMCMC(
            # graph_file, cooling_period=cooling_period, verbose=verbose)

        # Generate alternative plans
        if i < cooling_period:
            mcmc.generate_alternative_plans(rounds, rank, i)
        elif i >= cooling_period and i % 10 == 0:
            jjson = mcmc.generate_alternative_plans(rounds, rank, i)
            # mcmc.generate_alternative_plans(rounds)
            # Plot the data for the results of the recombinations
            # mcmc.plot_data()
            output.append(jjson)
        else:
            mcmc.generate_alternative_plans(rounds, rank, i)
    pdf = pd.DataFrame(output, columns=['plan'])
    fname = "update/recombination_of_districtsz" + str(rank)  # checkadd
    with open(fname, 'w') as file:
        file.write((pdf.to_json()))
    file.close()
    MD_refined = pd.read_json(
        './src/data/VA_refined_final.json')  # checkadd
    read_f = 'update/recombination_of_districtsz'+str(rank)
    df_final = pd.read_json(read_f)  # checkadd
    recomb = df_final
    precincts = []
    population = []
    white = []
    black = []
    asian = []
    hisp = []
    amin = []
    nhpi = []
    vap = []
    wvap = []
    bvap = []
    asianvap = []
    hvap = []
    aminvap = []
    nhpivap = []

    num_dist = 11  # total number of districts # 8 for md 11 for va
    print(len(recomb))
    for i in range(len(recomb)):  # total number of districting plans
        pprecincts = []
        ppopulation = []
        pwhite = []
        pblack = []
        pasian = []
        phisp = []
        pamin = []
        pnhpi = []
        pvap = []
        pwvap = []
        pbvap = []
        pasianvap = []
        phvap = []
        paminvap = []
        pnhpivap = []
        res2 = json.loads(recomb.plan[i])
        for j in range(num_dist):
            tprecincts = []
            tpopulation = 0
            twhite = 0
            tblack = 0
            tasian = 0
            thisp = 0
            tamin = 0
            tnhpi = 0
            tvap = 0
            twvap = 0
            tbvap = 0
            tasianvap = 0
            thvap = 0
            taminvap = 0
            tnhpivap = 0
            for k in range(len(res2.get("nodes"))):
                if int(res2.get("nodes")[k].get("district")) == j+1:
                    tprecincts.append(int(res2.get("nodes")[k].get("id")))
                    tpopulation = tpopulation + \
                        int(res2.get("nodes")[k].get("population"))
                    for l in range(2439):  # total number of precincts 2439 va 1809 md
                        if l == int(res2.get("nodes")[k].get("id")):
                            twhite = twhite + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("WHITE"))
                            tblack = tblack + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("BLACK"))
                            tasian = tasian + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("ASIAN"))
                            thisp = thisp + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("HISP"))
                            tamin = tamin + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("AMIN"))
                            tnhpi = tnhpi + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("NHPI"))
                            tvap = tvap + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("VAP"))
                            twvap = twvap + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("WVAP"))
                            tbvap = tbvap + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("BVAP"))
                            tasianvap = tasianvap + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("ASIANVAP"))
                            thvap = thvap + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("HVAP"))
                            taminvap = taminvap + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("AMINVAP"))
                            tnhpivap = tnhpivap + \
                                int(MD_refined.md_f.features[l].get(
                                    "properties").get("NHPIVAP"))
            pprecincts.append(tprecincts)
            ppopulation.append(tpopulation)
            pwhite.append(twhite)
            pblack.append(tblack)
            pasian.append(tasian)
            phisp.append(thisp)
            pamin.append(tamin)
            pnhpi.append(tnhpi)
            pvap.append(tvap)
            pwvap.append(twvap)
            pbvap.append(tbvap)
            pasianvap.append(tasianvap)
            phvap.append(thvap)
            paminvap.append(taminvap)
            pnhpivap.append(tnhpivap)
        precincts.append(pprecincts)
        population.append(ppopulation)
        white.append(pwhite)
        black.append(pblack)
        asian.append(pasian)
        hisp.append(phisp)
        amin.append(pamin)
        nhpi.append(pnhpi)
        vap.append(pvap)
        wvap.append(pwvap)
        bvap.append(pbvap)
        asianvap.append(pasianvap)
        hvap.append(phvap)
        aminvap.append(paminvap)
        nhpivap.append(pnhpivap)

    key = ["districtNumber", "population", "white", "black", "asian", "hisp", "amin",
           "nhpi", "vap", "wvap", "bvap", "asianvap", "hvap", "aminvap", "nhpivap", "precincts"]
    plans = []
    for j in range(len(population)):
        districts = []
        for i in range(len(population[0])):
            d1 = []
            d1.append(i+1)
            d1.append(population[j][i])
            d1.append(white[j][i])
            d1.append(black[j][i])
            d1.append(asian[j][i])
            d1.append(hisp[j][i])
            d1.append(amin[j][i])
            d1.append(nhpi[j][i])
            d1.append(vap[j][i])
            d1.append(wvap[j][i])
            d1.append(bvap[j][i])
            d1.append(asianvap[j][i])
            d1.append(hvap[j][i])
            d1.append(aminvap[j][i])
            d1.append(nhpivap[j][i])
            d1.append(precincts[j][i])
            it = iter(d1)
            dic1 = dict(zip(key, it))
            districts.append(dic1)
        dict2 = {"disctricts": districts}
        plans.append(dict2)
    dpf = {'plans': plans}
    dp = pd.DataFrame(dpf)
    output = []
    for i in range(len(dp)):
        output.append(dp['plans'][i])
    od = {'plans': output}
    fname = 'output/va-output'+str(rank)  # checkadd
    with open(fname, 'w') as f:
        json.dump(od, f)
    f.close()


# Save the data
if __name__ == "__main__":
    main()
