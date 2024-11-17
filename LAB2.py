
# pip install requests pandas networkx matplotlib streamlit

import requests
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import streamlit as st

# from BioGRID
def retrieve_ppi_biogrid(target_protein):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "66152c436cf41950ebdbf871e6ce284a",  # Replace with your BioGRID API key
        "format": "json",
        "searchNames": True,
        "geneList": target_protein,
        "organism": 9606,
        "searchbiogridids":True,
        "includeInteractors": True
    }
    response = requests.get(biogrid_url, params=params)
    if response.status_code != 200:
        st.error(f"Failed to retrieve data from BioGRID (Status Code: {response.status_code}).")
        return pd.DataFrame()  # Return empty DataFrame on failure
    
    network = response.json()
    
    if not network:  # Empty JSON response check
        st.warning(f"No interaction data found in BioGRID for protein: {target_protein}. Please check the protein name and try again.")
        return pd.DataFrame()
    network_df = pd.DataFrame.from_dict(network, orient='index')
    network_df.rename(columns={"OFFICIAL_SYMBOL_A": "Protein_A", "OFFICIAL_SYMBOL_B": "Protein_B"}, inplace=True)
    return network_df


# from STRING
def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606
    }

    response = requests.get(string_url, params=params)
    if response.status_code != 200:
        st.error("Failed to retrieve data from STRING.")
        return pd.DataFrame()  # Return empty DataFrame on failure
    network = response.json()
    network_df = pd.json_normalize(network)
    if not network_df.empty:
        # Rename columns for consistency
        network_df.rename(columns={"preferredName_A": "Protein_A", "preferredName_B": "Protein_B"}, inplace=True)
    return network_df

def generate_network(network_df):
    network_graph = nx.from_pandas_edgelist(network_df, 'Protein_A', 'Protein_B')
    return network_graph


def get_centralities(network_graph):
    centralities = {
        "Degree Centrality": nx.degree_centrality(network_graph),
        "Betweenness Centrality": nx.betweenness_centrality(network_graph),
        "Closeness Centrality": nx.closeness_centrality(network_graph),
        "Eigenvector Centrality": nx.eigenvector_centrality(network_graph),
        "PageRank Centrality": nx.pagerank(network_graph)
    }
    return centralities

# Streamlit app

# Title and description
st.title("Protein-Protein Interaction (PPI) Network Analysis")
st.write("This application allows you to retrieve Protein-Protein Interaction (PPI) networks and calculate centralities.")

# User input: Protein ID and selection of database
target_protein = st.text_input("Enter Target Protein(s)")
source = st.selectbox("Select PPI Database", ("BioGRID", "STRING"))

if target_protein:
    # Retrieve PPI data based on the selected source
    if source == "BioGRID":
        network_df = retrieve_ppi_biogrid(target_protein)
    elif source == "STRING":
        network_df = retrieve_ppi_string(target_protein)

    if not network_df.empty:
        # Generate the network graph
        network_graph = generate_network(network_df)

        # Get centrality measures
        centralities = get_centralities(network_graph)

        # Display results in two columns
        col1, col2 = st.columns(2)

        with col1:
            # Display PPI Data Information
            st.subheader("PPI Data Information")
            st.dataframe(network_df)  # Display the PPI DataFrame
            st.write(f"Number of Nodes: {network_graph.number_of_nodes()}")
            st.write(f"Number of Edges: {network_graph.number_of_edges()}")

            # Display network graph visualization
            st.subheader("Network Graph Visualization")
            plt.figure(figsize=(10, 10))
            nx.draw(network_graph, with_labels=True, node_color='skyblue', node_size=3000, font_size=10, font_weight='bold')
            st.pyplot(plt)
            plt.close()  # Close the plot to prevent overlap on reruns

        with col2:
            # Display Centrality Measures
            st.subheader("Centrality Measures")
            for centrality_name, centrality_value in centralities.items():
                st.write(f"{centrality_name}:")
                centrality_df = pd.DataFrame(list(centrality_value.items()), columns=["Node", centrality_name])
                st.dataframe(centrality_df)
    else:
        st.warning("No data available for the given protein.")