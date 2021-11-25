from pathlib import Path
import json
import math
import random
import networkx as nx
from networkx.algorithms import tree, boundary
import matplotlib.pyplot as plt
from functools import reduce
from networkx.readwrite import json_graph
import argparse

import pandas as pd
#from shapely.ops import unary_union
#from shapely.geometry import Polygon, MultiPolygon

# Globals we care about
MAX_POP_DIFFERENCE_PERCENTAGE = .30


class GerrymanderingMCMC():
    """
        Use a Markov Chain Monte Carlo simulation to generate an ensemble of redistricting plans
        against a given potential plan, and use the alternatives to see how much of an outlier the proposed plan is.
    """

    def __init__(self, graph_file, cooling_period=50, rounds=50, verbose=False):
        # We initialize all_districts here, but we really establish it when we read our graph in
        self.all_districts = set()
        self.g = self.read_graph(graph_file)
        self.cooling_period = cooling_period
        self.verbose = verbose
        self.district_colors = {
            "A": "red",
            "B": "green",
            "C": "yellow",
            "D": "blue"
        }
        self.data = []
        self.original_data = {}
        self.__record_key_stats(self.g, is_original_plan=True)

    def read_graph(self, path):
        """
            Given a path to an specialized JSON format describing a graph and it's metadata
            Returns a nx.Graph with relevant metadata stored on the node objects
        """
        g = nx.Graph()
        node_data = self.__load_json(path)
        for node_label in node_data:
            # Read in nodes from json format
            g.add_node(node_label)
            for adj_node in node_data[node_label]["adjacent_nodes"]:
                g.add_edge(node_label, adj_node)

            g.nodes[node_label]["population"] = node_data[node_label]["population"]

            #g.nodes[node_label]['VAP'] = node_data[node_label]['VAP']
            #g.nodes[node_label]['NH_WHITE'] = node_data[node_label]['NH_WHITE']
            #g.nodes[node_label]['WVAP'] = node_data[node_label]['WVAP']
            #g.nodes[node_label]['NH_BLACK'] = node_data[node_label]['NH_BLACK']
            #g.nodes[node_label]['BVAP'] = node_data[node_label]['BVAP']
            #g.nodes[node_label]['NH_ASIAN'] = node_data[node_label]['NH_ASIAN']
            #g.nodes[node_label]['ASIANVAP'] = node_data[node_label]['ASIANVAP']
            #g.nodes[node_label]['HISP'] = node_data[node_label]['HISP']
            #g.nodes[node_label]['HVAP'] = node_data[node_label]['HVAP']

            g.nodes[node_label]["voting_history"] = node_data[node_label]["voting_history"]
            #g.nodes[node_label]["geometry"] = node_data[node_label]["geometry"]
            g.nodes[node_label]["district"] = node_data[node_label]["district"]
            self.all_districts.add(node_data[node_label]["district"])

        fname = "output/original"
        with open(fname, 'w') as file:
            file.write(json.dumps(json_graph.adjacency_data(g)))

        return g

    def __load_json(self, path):
        """
            Loads a json file at a particular file-path location (str)
            Returns a json object corresponding to the data at the specified path
        """
        with Path(path).open() as json_file:
            return json.load(json_file)

    def __get_node_colors(self, g):
        """
            Use self.district_colors_map to return a list of colors corresponding to each node in a graph
        """
        return [self.district_colors[g.nodes[n]["district"]] for n in g.nodes]

    def __efficiency_gap(self, graph):
        """
            Determines the efficiency gap for a given district, and for whom it is in favor
            NOTE: Currently assumes a two-party system because plurality voting is the status quo; I'll improve the code when we improve our voting system
        """
        d_votes_wasted = r_votes_wasted = 0.00

        # Initialize tally-counters for both parties in each district
        district_dict = {}
        for district_label in self.all_districts:
            district_dict[district_label] = {"D": 0, "R": 0}

        # Count the votes for each party in each district
        for precinct_label in graph.nodes:
            precinct = graph.nodes[precinct_label]
            precinct_district = precinct["district"]
            # district_dict[precinct_district][precinct["voting_history"]] += 1

        # For each district, determine which party won and the wasted votes accordingly
        for district_label in self.all_districts:
            district = self.__get_district_subgraph(graph, district_label)
            # Total vote for a representative is just the number of precincts there are - the number of nodes on the graph
            total_district_votes = len(district.nodes)
            # Plurality voting would mean that we wouldn't even need this many votes if more than two systems were relevant players;
            #   but plurality voting also has a funny way of pushing elections towards two-party systems - so let's just assume that
            #   only two parties matter and therefore 50% of the precincts are needed
            votes_to_win = math.ceil(total_district_votes / 2.0)
            d_votes = district_dict[district_label]["D"]
            r_votes = district_dict[district_label]["R"]

            # And again - votes are wasted if they are cast for a losing party, or are a surplus beyond the amount required to win
            if d_votes > r_votes:
                d_votes_wasted += d_votes - votes_to_win
                r_votes_wasted += r_votes
            elif r_votes > d_votes:
                d_votes_wasted += d_votes
                r_votes_wasted += r_votes - votes_to_win
            else:
                None
                # NOTE: Pass when the district ends in a tie
        return (max([d_votes_wasted, r_votes_wasted]) - min([d_votes_wasted, r_votes_wasted])) / len(graph.nodes)

    def __random_district_label(self, graph):
        """
            UTIL method for getting a random district from a graph of nodes
        """
        return random.sample(graph, 1)[0]

    def __find_neighboring_district(self, g, district):
        """
            Given a graph and a district,
            Return a district_subgraph that neighbors it
        """
        node_boundary = boundary.node_boundary(g, district.nodes)
        districts_on_boundary = set(
            [g.nodes[n]["district"] for n in node_boundary])
        random_district_label = self.__random_district_label(
            districts_on_boundary)
        return self.__get_district_subgraph(g, random_district_label)

    def __get_district_nodes(self, g, district_label):
        """
            Given a nx.graph and a district_label
            Return a list of all the precincts (nodes) in that district
        """
        return [n for n in g.nodes if g.nodes[n]["district"] == district_label]

    def __get_district_subgraph(self, g, district_label):
        """
            Given a nx.graph and a district_label
            Return a subgraph view (not clone) corresponding to the precincts in that district
        """
        relevant_nodes = self.__get_district_nodes(g, district_label)
        return g.subgraph(relevant_nodes)

    def __district_size(self, potential_district):
        """
            Given a potential district of nodes, using the population size of the district
            Return the population size of the district
        """
        return reduce(lambda total, precinct: total + int(potential_district.nodes[precinct]["population"]),
                      potential_district.nodes(), 0)

    def __is_valid_district_plan(self, edge, mst_combined_subgraph, g):
        """
            For a given potential edge cut on an MST in the ReCom algorithm,
            Determine whether a series of required conditions is satisfied, including:
                1. Population size after new districting,
        """
        # 1. Check the population size after this cut
        # Does cutting this edges create two components with similar population sizes?
        (tail, head) = edge
        mst_combined_subgraph.remove_edge(tail, head)
        components = list(nx.connected_components(mst_combined_subgraph))
        comp_1 = g.subgraph(components[0])
        comp_2 = g.subgraph(components[1])
        pop_total = abs(self.__district_size(comp_1) +
                        self.__district_size(comp_2))
        pop_diff = abs(self.__district_size(comp_1) -
                       self.__district_size(comp_2))
        # Add edge back in case this doesn't work
        mst_combined_subgraph.add_edge(tail, head)
        return pop_diff < (MAX_POP_DIFFERENCE_PERCENTAGE * pop_total)

    def __update_new_districts_with_cut(self, edge, mst_combined_subgraph, g, d1, d2):
        """
            After chcecking that an edge chould be cut to create new districts after combining into a single mega-district
            Redistrict the new components accordingly
        """
        (tail, head) = edge
        mst_combined_subgraph.remove_edge(tail, head)
        components = list(nx.connected_components(mst_combined_subgraph))
        comp_1 = g.subgraph(components[0])
        comp_2 = g.subgraph(components[1])
        for node in comp_1.nodes:
            g.nodes[node]["district"] = d1
        for node in comp_2.nodes:
            g.nodes[node]["district"] = d2

    def recombination_of_districts(self, i):
        """
            Given a graph
            Perform the recombination algorithm described in https://mggg.org/va-report.pdf
            ("Comparison of Districting Plans for the Virginia House of Delfates"; Metric Geometry and Gerrymandering Group; Section 2.3.2)
            Alternative resource: the recombination algorithm described in https://arxiv.org/pdf/1911.05725.pdf
            ("Recombination: A family of Markov chains for redistricting"; Daryl DeFord, Moon Duchin, and Justin Solomon)
        """
        # print("hello")
        graph = self.g.copy()
        # Randomly sample a district
        d1_label = self.__random_district_label(self.all_districts)
        d1 = self.__get_district_subgraph(graph, d1_label)
        # Select one of its neighboring districts
        d2 = self.__find_neighboring_district(graph, d1)
        d2_label = d2.nodes[self.__random_district_label(d2.nodes)]["district"]
        combined_subgraph = graph.subgraph(list(d1.nodes) + list(d2.nodes))
        cuttable = False
        attempt_count = 0
        while cuttable is False:
            mst_combined_subgraph = self.__random_spanning_tree(
                combined_subgraph)
            # For all edges in the MST
            for edge in mst_combined_subgraph.edges:
                # If cutting this edge produces a valid districting
                cond = self.__is_valid_district_plan(
                    edge, mst_combined_subgraph, graph)
                if (cond):
                    cuttable = True
                    self.__update_new_districts_with_cut(
                        edge, mst_combined_subgraph, graph, d1_label, d2_label)
                    return graph
                if (attempt_count == 1000):
                    print(
                        "WARNING: Failed to make a recom after > 1000 iterations") if self.verbose else None
                    return graph
                else:
                    attempt_count += 1
                    (tail, head) = edge
                    mst_combined_subgraph.add_edge(tail, head)

    def __drawGraph(self, G, fname, options=None):
        """
            Given a graph, plot and draw it with these default options, or with options
        """

        node_colors = self.__get_node_colors(G)
        options = options if options != None else {
            'node_color': node_colors,
            'node_size': 100,
            "with_labels": True,
            'width': 3,
            "font_weight": "bold"
        }
        plt.figure(figsize=(8, 8))

        # Draw the graph we've been provided
        nx.draw_networkx(G, **options)
        # plt.show()
        fname = fname + ".png"
        plt.savefig(fname)
        plt.clf()

    def __random_spanning_tree(self, graph):
        """
            Given a graph
            Return a random spanning tree
        """
        for edge in graph.edges:
            graph.edges[edge]["weight"] = random.random()

        spanning_tree = tree.maximum_spanning_tree(
            graph, algorithm="kruskal", weight="weight"
        )
        return spanning_tree

    def __record_key_stats(self, graph, is_original_plan=False):
        """
            Given a potential districting plan (graph) and an optional flag for saying this is the original plan,
            Update our local data record to include stats for this plan
        """
        data_obj = {}
        # data_obj["eg"] = self.__efficiency_gap(graph)
        data_obj["d_districts"] = self.__count_votes(graph, "D")
        data_obj["r_districts"] = self.__count_votes(graph, "R")
        if is_original_plan:
            self.original_data = data_obj
        else:
            self.data.append(data_obj)

    def __winning_party_for_district(self, graph, district_label):
        """
            Given a graph and a district label,
            Return the party with the most precint votes
        """
        district = self.__get_district_subgraph(graph, district_label)
        # TODO: Update to use more than two parties
        demo_count = reduce(lambda demo_count, n_label: demo_count + 1 if district.nodes[n_label][
            "voting_history"] == "D" else demo_count - 1,
            district.nodes, 0)
        if demo_count == 0:
            return None
        elif demo_count < 0:
            return "R"
        else:
            return "D"

    def __count_votes(self, graph, party):
        """
            Given a graph and party,
            Return the number of districts that voted for that party
        """
        return reduce(
            lambda count, d_label: count +
            1 if self.__winning_party_for_district(
                graph, d_label) == party else count,
            self.all_districts, 0)

    def plot_data(self):
        print("")
        # plt.hist([d["eg"] for d in self.data], bins="auto", alpha=0.5, facecolor='blue')
        # plt.title("Efficiency Gap")
        # plt.axvline(self.original_data["eg"], label="Original Plan")
        # plt.legend()
        # # plt.show()
        # plt.savefig('efficiency_gap.png')
        # plt.clf()
        # plt.hist([d["d_districts"] for d in self.data], bins=5, range=(0,4), alpha=0.5, facecolor='blue')
        # plt.title("Democratic Districts")
        # plt.axvline(self.original_data["d_districts"], label="Original Plan")
        # plt.legend()
        # # plt.show()
        # plt.savefig('democratic_districts.png')
        # plt.clf()
        # plt.hist([d["r_districts"] for d in self.data], bins=5, range=(0,4), alpha=0.5, facecolor='blue')
        # plt.title("Republican Districts")
        # plt.axvline(self.original_data["r_districts"], label="Original Plan")
        # plt.legend()
        # # plt.show()
        # plt.savefig('republican_districts.png')
        # plt.clf()

    def generate_district_geometry(self, graph):
        # Add district_geometries field to the JSON object
        graph["district_geometries"] = []
        # Retrieve precincts from graph JSON Object
        precincts_df = pd.json_normalize(graph, record_path=['nodes'])
        district_groups = precincts_df.groupby(by='district')

        # For each district group, form a Polygon union.
        for district_name, district_df in district_groups:
            # Initialize a list to save precinct polygons
            polygons = list()

            # Iterate through each precinct and save its coordinates.
            for index, row in district_df.iterrows():
                if row['geometry.type'] == 'Polygon':
                    polygons.append(Polygon(row['geometry.coordinates'][0]))
                else:
                    for coords in row['geometry.coordinates']:
                        polygons.append(Polygon(coords[0]))

            # Join these precincts coordinates to form district geometry
            union = unary_union(polygons)

            # Construct district with appropriate key-value pairs.
            if isinstance(union, MultiPolygon):
                geom_type = 'MultiPolygon'
                coordinates = [[list(poly.exterior.coords)]
                               for poly in list(union)]
            else:
                geom_type = 'Polygon'
                coordinates = [list(union.exterior.coords)]

            # Update district_geometries field of passed graph JSON Object
            graph["district_geometries"].append(dict(district_name=district_name,
                                                     totalPop=int(
                                                         district_df['population'].sum()),
                                                     totalVAP=int(
                                                         district_df['VAP'].sum()),
                                                     totalNH_WHITE=int(
                                                         district_df['NH_WHITE'].sum()),
                                                     totalWVAP=int(
                                                         district_df['WVAP'].sum()),
                                                     totalNH_BLACK=int(
                                                         district_df['NH_BLACK'].sum()),
                                                     totalBVAP=int(
                                                         district_df['BVAP'].sum()),
                                                     totalNH_ASIAN=int(
                                                         district_df['NH_ASIAN'].sum()),
                                                     totalASIANVAP=int(
                                                         district_df['ASIANVAP'].sum()),
                                                     totalHISP=int(
                                                         district_df['HISP'].sum()),
                                                     totalHVAP=int(
                                                         district_df['HVAP'].sum()),
                                                     geometry=dict(type=geom_type, coordinates=coordinates)))
        return graph

    def generate_recur_input(self, rank):
        # print('hi')
        ref = "update/update"+str(rank)
        a = pd.read_json(ref)  # checkadd
        a_load = json.loads(a.plan[0])
        aa = []
        bb = []
        for i in range(len(a_load.get('nodes'))):
            aa.append(0)
            bb.append(0)
        for i in range(len(aa)):
            ind = int(a_load.get('nodes')[i].get('id'))
            aa[ind] = a_load.get('nodes')[i]
            bb[ind] = a_load.get('adjacency')[i]
        neighbors = []
        for i in range(len(bb[0])):
            neighbors.append(int(bb[0][i].get('id')))
        neighbors = sorted(neighbors)
        res = []
        for j in range(len(neighbors)):
            res.append("% s" % neighbors[j])
        data = ['0', aa[0].get('population'), "D", aa[0].get('district')]
        col = ['adjacent_nodes', 'population', 'voting_history', 'district']
        mg_df = pd.DataFrame(data, col)
        mg_df.at['adjacent_nodes', 0] = res
        for i in range(len(a_load.get('nodes'))-1):
            data = ['0', aa[i+1].get('population'), "D",
                    aa[i+1].get('district')]
            mg_df[i+1] = data
            neighbors = []
            # print(i)
            for j in range(len(bb[i+1])):
                if int(bb[i+1][j].get('id')) > i+1:
                    neighbors.append(int(bb[i+1][j].get('id')))
            neighbors = sorted(neighbors)
            res = []
            for j in range(len(neighbors)):
                res.append("% s" % neighbors[j])
            # print(res)
            mg_df.at['adjacent_nodes', i+1] = res
        jjson = pd.DataFrame.to_json(mg_df)
        fname_update = "./src/data/prep2_va"+str(rank)  # checkadd
        with open(fname_update, 'w') as fo:
            fo.write(jjson)
        fo.close()

    def generate_alternative_plans(self, rounds, rank, ind):
        print(ind)
        # Run `cooling`-many rounds to randomize the plan a bit
        outputf = []
        fname_update = "update/update"+str(rank)  # checkadd
        #fname_update = "output/update"
        # self.__drawGraph(self.g, "output/original")
        if ind < self.cooling_period:
            # save to update and read
            print("Randomizing the seed plan",
                  ind) if ind % 25 == 0 and self.verbose else None
            graph = self.recombination_of_districts(ind)
            outputf.append(json.dumps(json_graph.adjacency_data(graph)))
            # self.__drawGraph(graph, fname)
            pdf = pd.DataFrame(
                outputf, columns=['plan'])
            # output
            with open(fname_update, 'w') as file:
                file.write((pdf.to_json()))
            file.close()
            # update (output to input)
            self.generate_recur_input(rank)

        else:
            # Run `rounds`-many recombinations to build a distribution of a few key stats
            # save to update and read and return json
            print("Finding recomb ... ",
                  ind) if ind % 20 == 0 and self.verbose else None
            graph = self.recombination_of_districts(ind)
            outputf.append(json.dumps(json_graph.adjacency_data(graph)))

            self.__record_key_stats(graph)
            #fname = "output/recombination_of_districts" + str(rank)
            # self.__drawGraph(graph, fname)
            pdf = pd.DataFrame(
                outputf, columns=['plan'])
            # output
            with open(fname_update, 'w') as file:
                file.write((pdf.to_json()))
            file.close()
            self.generate_recur_input(rank)
            return json.dumps(json_graph.adjacency_data(graph))

        print("DONE 123 Finding alternative district plans") if self.verbose else None
