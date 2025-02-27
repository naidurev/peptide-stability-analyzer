# Essential imports that are available in Pyodide
import numpy as np
from typing import List, Tuple, Dict
import json

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
        if cys_count == 0:
            raise ValueError("No cysteines found in the sequence")

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

def evaluate_pattern_stability(self, pattern: List[Tuple[int, int]]) -> Dict:
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

def _get_charge_positions(self) -> List[Tuple[int, float]]:
        """Get positions and charges of ionizable residues"""
        charges = []
        charge_map = {'K': 1, 'R': 1, 'H': 0.5, 'D': -1, 'E': -1}
        for i, aa in enumerate(self.sequence):
            if aa in charge_map:
                charges.append((i, charge_map[aa]))
        return charges
		
def _calculate_solvation_energy(self, pattern: List[Tuple[int, int]]) -> float:
        """Calculate solvation energy using empirical transfer free energies"""
        energy = 0.0
        burial_state = self._estimate_burial_state(pattern)

        for i, aa in enumerate(self.sequence):
            if aa in self.energy_params['solvation']:
                energy += self.energy_params['solvation'][aa] * burial_state[i]
        return energy

def _estimate_burial_state(self, pattern: List[Tuple[int, int]]) -> List[float]:
        """Estimate burial state (0-1) for each residue"""
        burial = [0.0] * self.length
        for cys1, cys2 in pattern:
            for i in range(cys1, cys2+1):
                burial[i] = max(burial[i], 0.5)
            burial[cys1] = burial[cys2] = 0.8
        return burial

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
            loop_size = abs(cys2 - cys1)
            # Simplified calculation based on loop size
            if loop_size > 0:
                sigma = 4.0
                epsilon = 0.2
                # Approximation of average distance
                distance = max(3.0, loop_size * 0.5)
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