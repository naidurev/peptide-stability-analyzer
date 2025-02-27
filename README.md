# Peptide Stability Analyzer

A web-based tool for analyzing peptide stability and predicting optimal disulfide bond patterns. This application helps researchers and scientists evaluate different disulfide bonding configurations and their impact on peptide stability.

## Features

- Analyze peptide sequences with cysteines to identify possible disulfide bonding patterns
- Calculate energy components that contribute to stability:
  - Conformational energy
  - Electrostatic interactions
  - Solvation energy
  - Strain energy
  - Entropy
  - Van der Waals forces
- Visualize results with interactive charts and tables
- Works entirely in the browser with no server-side processing
- Responsive design that works on desktop and mobile devices

## How to Use

1. Visit the application in your web browser
2. Enter your peptide sequence using standard one-letter amino acid codes
   - Example: `CGKCTLKECHPSKGCSR`
   - Your sequence must contain an even number of cysteines (C)
3. Click "Analyze Stability" to generate results
4. View results in three tabs:
   - **Energy Components**: Chart showing energetic contributions for top patterns
   - **Pattern Analysis**: Table of all possible disulfide patterns ranked by stability
   - **Structure Visualization**: Visual representation of your sequence with best pattern

## Technical Details

This application uses:
- Pyodide to run Python scientific code directly in the browser
- Chart.js for data visualization
- Bootstrap for responsive design
- NumPy for scientific calculations

The peptide stability calculations are based on:
- Conformational preferences of disulfide bonds
- Electrostatic interactions using Debye-Hückel theory
- Solvation energies based on amino acid transfer free energies
- Strain energy from bond geometry
- Entropic effects from loop formation

## Setup Instructions

### Option 1: Run directly from GitHub Pages
The application is hosted on GitHub Pages and can be accessed directly at [URL].

### Option 2: Run locally
1. Clone this repository: `git clone https://github.com/yourusername/peptide-stability-analyzer.git`
2. Navigate to the project folder: `cd peptide-stability-analyzer`
3. Serve the files using a local web server. For example:
   - Using Python: `python -m http.server`
   - Using Node.js: `npx serve`
4. Open your browser and navigate to `http://localhost:8000` (or the port specified by your server)

## Requirements

- Modern web browser with JavaScript enabled (Chrome, Firefox, Edge, or Safari)
- No installation required - everything runs in the browser

## Troubleshooting

- **Loading fails**: Ensure you have a stable internet connection as the application needs to download the Pyodide environment
- **Analysis takes too long**: Try with a shorter peptide sequence, as computation time increases with sequence length
- **Browser crashes**: If using a very large sequence, try using a desktop browser with more memory

## License

MIT License
