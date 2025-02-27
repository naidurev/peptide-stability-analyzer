import numpy as np
from Bio import SeqUtils
from Bio.PDB import *
from typing import List, Tuple, Dict
import ipywidgets as widgets
from IPython.display import display, HTML
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

class PeptideStabilityAnalyzer:
    def __init__(self, sequence: str):
        """
        Initialize analyzer with peptide sequence

        Args:
            sequence: Amino acid sequence in single letter code
        """
        self.sequence = sequence
        self.length = len(sequence)
        self.cys_positions = [i for i, aa in enumerate(sequence) if aa == 'C']

        # Energy parameters
        self.energy_params = {
            'torsional': {
                'chi1': {'preferred': [-60, 180], 'penalty': 2.0},
                'chi2': {'preferred': [-60, 60], 'penalty': 1.5},
                'chi3': {'preferred': [-120, 120], 'penalty': 1.0}
            },
            'electrostatic': {
                'dielectric': 80.0,  # Water dielectric constant
                'screening_length': 10.0  # Debye screening length (Å)
            },
            'solvation': {
                # Transfer free energies (kcal/mol)
                'ASP': -10.95, 'GLU': -10.20, 'HIS': -4.66, 'LYS': -9.52,
                'ARG': -10.28, 'CYS': -1.24, 'MET': 2.35, 'TYR': -0.14,
                'TRP': 3.07, 'PHE': 3.79, 'ALA': 1.94, 'ILE': 4.92,
                'LEU': 4.92, 'PRO': -0.99, 'VAL': 4.04, 'GLY': 0.94,
                'SER': -5.06, 'THR': -4.88, 'ASN': -9.68, 'GLN': -9.38
            }
        }

    def get_possible_disulfide_patterns(self) -> List[List[Tuple[int, int]]]:
        """Generate all possible disulfide bonding patterns"""
        cys_count = len(self.cys_positions)
        if cys_count % 2 != 0:
            raise ValueError("Odd number of cysteines cannot form complete disulfide patterns")

        patterns = []
        def generate_patterns(remaining_pos, current_pattern):
            if not remaining_pos:
                patterns.append(current_pattern[:])
                return
            current = remaining_pos[0]
            for partner in remaining_pos[1:]:
                new_pattern = current_pattern + [(current, partner)]
                new_remaining = [p for p in remaining_pos if p != current and p != partner]
                generate_patterns(new_remaining, new_pattern)

        generate_patterns(self.cys_positions, [])
        return patterns

    def evaluate_pattern_stability(self, pattern: List[Tuple[int, int]], structure=None) -> Dict:
        """
        Evaluate stability contribution of a specific disulfide pattern with detailed energy terms
        """
        energy_components = {
            'conformational_energy': self._calculate_conformational_energy(pattern),
            'electrostatic_energy': self._calculate_electrostatic_energy(pattern),
            'solvation_energy': self._calculate_solvation_energy(pattern),
            'strain_energy': self._calculate_detailed_strain_energy(pattern),
            'entropy': self._calculate_detailed_entropy(pattern),
            'van_der_waals': self._calculate_vdw_energy(pattern)
        }

        total_energy = sum(energy_components.values())

        stability_metrics = {
            'energy_components': energy_components,
            'total_energy': total_energy,
            'loop_sizes': self._calculate_loop_sizes(pattern),
            'buried_surface_area': self._estimate_buried_surface(pattern)
        }

        return stability_metrics

    def _calculate_loop_sizes(self, pattern: List[Tuple[int, int]]) -> List[int]:
        """Calculate sizes of loops formed by disulfide bonds"""
        return [abs(cys2 - cys1) for cys1, cys2 in pattern]

    def _calculate_conformational_energy(self, pattern: List[Tuple[int, int]]) -> float:
        """Calculate conformational energy considering backbone and side-chain conformations"""
        energy = 0.0
        for cys1, cys2 in pattern:
            loop_size = abs(cys2 - cys1)
            if loop_size < 3:
                energy += 10.0
            elif loop_size > 15:
                energy += 0.5 * (loop_size - 15)

            # Side-chain conformational energy
            for chi, params in self.energy_params['torsional'].items():
                estimated_chi = self._estimate_chi_angle(cys1, cys2)
                if not (params['preferred'][0] <= estimated_chi <= params['preferred'][1]):
                    energy += params['penalty']
        return energy

    def _calculate_electrostatic_energy(self, pattern: List[Tuple[int, int]]) -> float:
        """Calculate electrostatic energy using Debye-Hückel theory"""
        energy = 0.0
        charges = self._get_charge_positions()

        for i, (pos1, q1) in enumerate(charges):
            for pos2, q2 in charges[i+1:]:
                distance = abs(pos2 - pos1) * 3.8
                if distance > 0:
                    energy += (q1 * q2 / (self.energy_params['electrostatic']['dielectric'] * distance) *
                             np.exp(-distance / self.energy_params['electrostatic']['screening_length']))
        return energy * 332.0

    def _calculate_solvation_energy(self, pattern: List[Tuple[int, int]]) -> float:
        """Calculate solvation energy using empirical transfer free energies"""
        energy = 0.0
        burial_state = self._estimate_burial_state(pattern)

        for i, aa in enumerate(self.sequence):
            if aa in self.energy_params['solvation']:
                energy += self.energy_params['solvation'][aa] * burial_state[i]
        return energy

    def _calculate_detailed_strain_energy(self, pattern: List[Tuple[int, int]]) -> float:
        """Calculate detailed strain energy including bond, angle, and torsional terms"""
        energy = 0.0
        for cys1, cys2 in pattern:
            ideal_bond_length = 2.05
            estimated_length = self._estimate_bond_length(cys1, cys2)
            energy += 100 * (estimated_length - ideal_bond_length)**2

            ideal_angle = 103
            estimated_angle = self._estimate_bond_angle(cys1, cys2)
            energy += 50 * (estimated_angle - ideal_angle)**2

            estimated_torsion = self._estimate_torsion_angle(cys1, cys2)
            energy += 10 * (1 + np.cos(3 * estimated_torsion))
        return energy

    def _calculate_vdw_energy(self, pattern: List[Tuple[int, int]]) -> float:
        """Calculate van der Waals energy using Lennard-Jones potential"""
        energy = 0.0
        for cys1, cys2 in pattern:
            loop_residues = self.sequence[cys1:cys2+1]
            for i, aa1 in enumerate(loop_residues):
                for j, aa2 in enumerate(loop_residues[i+1:], i+1):
                    distance = (j - i) * 3.8
                    if distance > 0:
                        sigma = 4.0
                        epsilon = 0.2
                        r6 = (sigma/distance)**6
                        energy += epsilon * (r6**2 - 2*r6)
        return energy

    def _calculate_detailed_entropy(self, pattern: List[Tuple[int, int]]) -> float:
        """Calculate detailed entropy changes"""
        entropy = 0.0
        for cys1, cys2 in pattern:
            loop_size = abs(cys2 - cys1)
            entropy -= 0.5 * loop_size
            entropy -= 2.0
            buried_residues = self._count_buried_residues(cys1, cys2)
            entropy += 0.2 * buried_residues
        return entropy

    def _get_charge_positions(self) -> List[Tuple[int, float]]:
        """Get positions and charges of ionizable residues"""
        charges = []
        charge_map = {'K': 1, 'R': 1, 'H': 0.5, 'D': -1, 'E': -1}
        for i, aa in enumerate(self.sequence):
            if aa in charge_map:
                charges.append((i, charge_map[aa]))
        return charges

    def _estimate_burial_state(self, pattern: List[Tuple[int, int]]) -> List[float]:
        """Estimate burial state (0-1) for each residue"""
        burial = [0.0] * self.length
        for cys1, cys2 in pattern:
            for i in range(cys1, cys2+1):
                burial[i] = max(burial[i], 0.5)
            burial[cys1] = burial[cys2] = 0.8
        return burial

    def _estimate_buried_surface(self, pattern: List[Tuple[int, int]]) -> float:
        """Estimate buried surface area"""
        return sum(abs(cys2 - cys1) for cys1, cys2 in pattern)

    def _estimate_chi_angle(self, cys1: int, cys2: int) -> float:
        """Estimate chi angle based on sequence context"""
        return -60 + (abs(cys2 - cys1) % 3) * 60

    def _estimate_bond_length(self, cys1: int, cys2: int) -> float:
        """Estimate disulfide bond length"""
        return 2.05

    def _estimate_bond_angle(self, cys1: int, cys2: int) -> float:
        """Estimate C-S-S bond angle"""
        return 103 + (abs(cys2 - cys1) % 2) * 2

    def _estimate_torsion_angle(self, cys1: int, cys2: int) -> float:
        """Estimate disulfide torsion angle"""
        return np.pi / 3

    def _count_buried_residues(self, cys1: int, cys2: int) -> int:
        """Count number of buried residues in a disulfide-bonded loop"""
        return len([i for i in range(cys1, cys2+1)
                   if self.sequence[i] in 'FILMVWY'])
def create_peptide_gui():
    """Create the GUI for peptide analysis"""
    style = {'description_width': 'initial'}

    sequence_input = widgets.Text(
        value='',
        placeholder='Enter amino acid sequence',
        description='Peptide Sequence:',
        style=style,
        layout={'width': '500px'}
    )

    analyze_button = widgets.Button(
        description='Analyze Stability',
        button_style='success',
        style=style
    )

    output_area = widgets.Output()

    tab_names = ['Energy Components', 'Pattern Analysis', 'Structure Visualization']
    tab = widgets.Tab([widgets.Output() for _ in tab_names])
    for i, name in enumerate(tab_names):
        tab.set_title(i, name)

    def on_analyze_button_clicked(b):
        sequence = sequence_input.value.upper()

        # Validate sequence
        if not sequence:
            with output_area:
                print("Please enter a peptide sequence.")
                return

        if not all(aa in 'ACDEFGHIKLMNPQRSTVWY' for aa in sequence):
            with output_area:
                print("Invalid amino acid sequence. Please use standard one-letter codes.")
                return

        # Clear previous outputs
        for tab_output in tab.children:
            tab_output.clear_output()

        try:
            analyzer = PeptideStabilityAnalyzer(sequence)
            patterns = analyzer.get_possible_disulfide_patterns()

            if not patterns:
                with output_area:
                    print("No valid disulfide patterns found (need even number of cysteines).")
                return

            # Analyze patterns
            results = []
            for pattern in patterns:
                stability_metrics = analyzer.evaluate_pattern_stability(pattern)
                results.append({
                    'pattern': pattern,
                    'metrics': stability_metrics
                })

            results.sort(key=lambda x: x['metrics']['total_energy'])

            # Tab 1: Energy Components
            with tab.children[0]:
                plt.figure(figsize=(12, 6))
                energy_data = []
                for i, result in enumerate(results[:5]):
                    components = result['metrics']['energy_components']
                    for component, value in components.items():
                        energy_data.append({
                            'Pattern': f'Pattern {i+1}',
                            'Component': component,
                            'Energy': value
                        })

                # Convert to DataFrame
                df_energy = pd.DataFrame(energy_data)

                # Create the plot
                sns.barplot(data=df_energy, x='Pattern', y='Energy', hue='Component')
                plt.title('Energy Components of Top 5 Patterns')
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.show()

            # Tab 2: Pattern Analysis
            with tab.children[1]:
                html_table = "<table style='width:100%; border-collapse: collapse;'>"
                html_table += "<tr><th>Pattern</th><th>Total Energy</th><th>Loop Sizes</th></tr>"

                for i, result in enumerate(results):
                    pattern_str = ' & '.join([f'C{c1+1}-C{c2+1}' for c1, c2 in result['pattern']])
                    html_table += f"<tr><td>{pattern_str}</td>"
                    html_table += f"<td>{result['metrics']['total_energy']:.2f}</td>"
                    html_table += f"<td>{', '.join(map(str, result['metrics']['loop_sizes']))}</td></tr>"

                html_table += "</table>"
                display(HTML(html_table))

            # Tab 3: Structure Visualization
            with tab.children[2]:
                plt.figure(figsize=(12, 3))
                sequence_list = list(sequence)
                colors = ['red' if aa == 'C' else 'black' for aa in sequence_list]

                plt.plot(range(len(sequence)), [0]*len(sequence), 'k-', alpha=0.2)
                for i, (aa, color) in enumerate(zip(sequence_list, colors)):
                    plt.text(i, 0, aa, color=color, ha='center', va='center', fontsize=12)

                # Plot best disulfide pattern
                best_pattern = results[0]['pattern']
                for cys1, cys2 in best_pattern:
                    plt.plot([cys1, cys2], [0.2, 0.2], 'b-', alpha=0.5)

                plt.title('Sequence Visualization with Best Disulfide Pattern')
                plt.axis('off')
                plt.tight_layout()
                plt.show()

        except Exception as e:
            with output_area:
                print(f"Error during analysis: {str(e)}")

    analyze_button.on_click(on_analyze_button_clicked)

    # Arrange widgets
    gui = widgets.VBox([
        widgets.HBox([sequence_input, analyze_button]),
        output_area,
        tab
    ])

    return gui

# Create and display the GUI
display(create_peptide_gui())

# Print usage instructions
print("\nPeptide Stability Analyzer Instructions:")
print("----------------------------------------")
print("1. Enter your peptide sequence using standard one-letter amino acid codes:")
print("   - Example sequence: CGKCTLKECHPSKGCSR")
print("   - Must contain an even number of cysteines (C)")
print("\n2. Click 'Analyze Stability' to perform the analysis")
print("\n3. View results in three tabs:")
print("   - Energy Components: Bar plot of energy terms for top patterns")
print("   - Pattern Analysis: Table of all possible disulfide patterns")
print("   - Structure Visualization: Sequence diagram showing best pattern")
print("\nValid amino acids: A C D E F G H I K L M N P Q R S T V W Y")
print("\nNote: Cysteines (C) are highlighted in red in the visualization")