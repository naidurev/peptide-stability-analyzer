// Global variables
let pyodide = null;
let analyzer = null;

// Initialize Pyodide when page loads
async function main() {
    try {
        // Show loading indicator
        document.getElementById('loading-indicator').style.display = 'flex';
        document.getElementById('app-container').style.display = 'none';
        
        // Load Pyodide
        pyodide = await loadPyodide();
        console.log("Pyodide loaded successfully");
        
        // Install required packages
        await pyodide.loadPackagesFromImports(`
            import numpy
            import pandas
            import matplotlib
        `);
        
        // Define our Python code (simplified version of the analyzer)
        const pythonCode = await fetchPythonCode();
        
        // Run the Python code to define our classes
        await pyodide.runPythonAsync(pythonCode);
        
        console.log("Python environment ready");
        
        // Hide loading indicator and show app
        document.getElementById('loading-indicator').style.display = 'none';
        document.getElementById('app-container').style.display = 'block';
        
        // Set up event listeners
        setupEventListeners();
        
    } catch (error) {
        console.error("Error initializing application:", error);
        document.getElementById('loading-indicator').innerHTML = `
            <div class="alert alert-danger">
                <h4>Error Loading Application</h4>
                <p>${error.message}</p>
                <p>Please try refreshing the page or using a different browser.</p>
            </div>
        `;
    }
}
// Fetch the Python code from the python file
async function fetchPythonCode() {
    try {
        const response = await fetch('python/peptide_analyzer.py');
        if (!response.ok) {
            throw new Error(`Failed to fetch Python code: ${response.status}`);
        }
        return await response.text();
    } catch (error) {
        console.error("Error fetching Python code:", error);
        throw error;
    }
}

// Set up event listeners for UI elements
function setupEventListeners() {
    const analyzeButton = document.getElementById('analyze-button');
    const sequenceInput = document.getElementById('sequence-input');
    
    // Handle analyze button click
    analyzeButton.addEventListener('click', () => {
        analyzeSequence(sequenceInput.value);
    });
    
    // Handle enter key in input field
    sequenceInput.addEventListener('keypress', (event) => {
        if (event.key === 'Enter') {
            analyzeSequence(sequenceInput.value);
        }
    });
    
    // Example sequence button (for demonstration)
    const exampleButton = document.createElement('button');
    exampleButton.textContent = 'Use Example';
    exampleButton.className = 'btn btn-outline-secondary ms-2';
    exampleButton.addEventListener('click', () => {
        sequenceInput.value = 'CGKCTLKECHPSKGCSR';
        analyzeSequence(sequenceInput.value);
    });
    
    // Add example button after the input group
    const inputGroup = sequenceInput.parentElement;
    inputGroup.insertAdjacentElement('afterend', exampleButton);
}
// Analyze the provided sequence
async function analyzeSequence(sequence) {
    // Validate sequence
    sequence = sequence.trim().toUpperCase();
    if (!sequence) {
        showStatusMessage('Please enter a peptide sequence.', 'danger');
        return;
    }
    
    const validAminoAcids = 'ACDEFGHIKLMNPQRSTVWY';
    if (![...sequence].every(aa => validAminoAcids.includes(aa))) {
        showStatusMessage('Invalid amino acid sequence. Please use standard one-letter codes.', 'danger');
        return;
    }
    
    try {
        // Call Python code to analyze the sequence
        const result = await pyodide.runPythonAsync(`
            try:
                # Create analyzer instance
                analyzer = PeptideStabilityAnalyzer("${sequence}")
                
                # Get possible patterns
                patterns = analyzer.get_possible_disulfide_patterns()
                
                # If no patterns, return empty result
                if not patterns:
                    import json
                    json.dumps({"error": "No valid disulfide patterns found (need even number of cysteines)."})
                else:
                    # Analyze patterns
                    results = []
                    for pattern in patterns:
                        stability_metrics = analyzer.evaluate_pattern_stability(pattern)
                        results.append({
                            "pattern": pattern,
                            "metrics": stability_metrics
                        })
                    
                    # Sort by total energy
                    results.sort(key=lambda x: x["metrics"]["total_energy"])
                    
                    # Convert to JSON
                    import json
                    json.dumps(results)
            except Exception as e:
                import json
                json.dumps({"error": str(e)})
        `);
		// Parse the JSON result
        const data = JSON.parse(result);
        
        // Check for errors
        if (data.error) {
            showStatusMessage(data.error, 'danger');
            return;
        }
        
        // Display results
        displayResults(data, sequence);
        showStatusMessage('Analysis complete!', 'success');
        
    } catch (error) {
        console.error("Error analyzing sequence:", error);
        showStatusMessage(`Error analyzing sequence: ${error.message}`, 'danger');
    }
}

// Display status message to the user
function showStatusMessage(message, type = 'info') {
    const statusMessage = document.getElementById('status-message');
    statusMessage.textContent = message;
    statusMessage.className = `alert alert-${type}`;
    statusMessage.style.display = 'block';
    
    // Hide message after 5 seconds if it's a success message
    if (type === 'success') {
        setTimeout(() => {
            statusMessage.style.display = 'none';
        }, 5000);
    }
}

// Display analysis results in the UI
function displayResults(results, sequence) {
    // Display top 5 patterns in energy chart
    createEnergyChart(results.slice(0, 5));
    
    // Display pattern table
    createPatternTable(results);
    
    // Display sequence visualization with best pattern
    createSequenceVisualization(sequence, results[0].pattern);
}
// Create energy components chart
function createEnergyChart(topResults) {
    const ctx = document.getElementById('energy-chart').getContext('2d');
    
    // Prepare data for chart
    const labels = topResults.map((_, index) => `Pattern ${index + 1}`);
    const datasets = [];
    
    // Get all energy component keys from the first result
    const componentKeys = Object.keys(topResults[0].metrics.energy_components);
    
    // Assign colors to each component
    const colors = [
        'rgba(54, 162, 235, 0.7)',  // blue
        'rgba(255, 99, 132, 0.7)',  // red
        'rgba(75, 192, 192, 0.7)',  // green
        'rgba(255, 206, 86, 0.7)',  // yellow
        'rgba(153, 102, 255, 0.7)', // purple
        'rgba(255, 159, 64, 0.7)'   // orange
    ];
    
    // Create a dataset for each energy component
    componentKeys.forEach((component, index) => {
        const data = topResults.map(result => result.metrics.energy_components[component]);
        
        datasets.push({
            label: component.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase()),
            data: data,
            backgroundColor: colors[index % colors.length],
            borderColor: colors[index % colors.length].replace('0.7', '1'),
            borderWidth: 1
        });
    });
    
    // Create and render the chart
    if (window.energyChart) {
        window.energyChart.destroy();
    }
    
    window.energyChart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: labels,
            datasets: datasets
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                x: {
                    stacked: true
                },
                y: {
                    stacked: true,
                    title: {
                        display: true,
                        text: 'Energy (arbitrary units)'
                    }
                }
            }
        }
    });
}
// Create pattern analysis table
function createPatternTable(results) {
    const tableBody = document.getElementById('pattern-table-body');
    tableBody.innerHTML = '';
    
    results.forEach((result, index) => {
        const row = document.createElement('tr');
        
        // Pattern column
        const patternCell = document.createElement('td');
        const patternStr = result.pattern.map(([c1, c2]) => `C${c1 + 1}-C${c2 + 1}`).join(' & ');
        patternCell.textContent = patternStr;
        row.appendChild(patternCell);
        
        // Total energy column
        const energyCell = document.createElement('td');
        energyCell.textContent = result.metrics.total_energy.toFixed(2);
        row.appendChild(energyCell);
        
        // Loop sizes column
        const loopCell = document.createElement('td');
        loopCell.textContent = result.metrics.loop_sizes.join(', ');
        row.appendChild(loopCell);
        
        // Highlight the best pattern
        if (index === 0) {
            row.classList.add('table-success');
        }
        
        tableBody.appendChild(row);
    });
}

// Create sequence visualization with disulfide bonds
function createSequenceVisualization(sequence, bestPattern) {
    const sequenceDisplay = document.getElementById('sequence-display');
    sequenceDisplay.innerHTML = '';
    
    // Display sequence with highlighted cysteines
    [...sequence].forEach((aa, index) => {
        const span = document.createElement('span');
        span.textContent = aa;
        span.classList.add('amino-acid');
        span.dataset.index = index;
        
        if (aa === 'C') {
            span.classList.add('cysteine');
        }
        
        sequenceDisplay.appendChild(span);
    });
    
    // Create the structure chart
    createStructureChart(sequence, bestPattern);
}
// Create structure visualization chart
function createStructureChart(sequence, bestPattern) {
    const ctx = document.getElementById('structure-chart').getContext('2d');
    
    // Convert best pattern to pairs of indices for visualization
    const bondPairs = bestPattern.map(([c1, c2]) => ({ from: c1, to: c2 }));
    
    // Create data points for the amino acids (on a straight line)
    const aaPositions = [...sequence].map((_, index) => ({
        x: index,
        y: 0
    }));
    
    // Create data for the chart
    const data = {
        labels: [...sequence].map((aa, idx) => `${aa}${idx + 1}`),
        datasets: [
            {
                type: 'scatter',
                label: 'Amino Acids',
                data: aaPositions,
                backgroundColor: [...sequence].map(aa => aa === 'C' ? 'red' : 'black'),
                pointRadius: 8,
                pointHoverRadius: 10
            }
        ]
    };
    
    // Add lines for disulfide bonds
    bondPairs.forEach((bond, index) => {
        data.datasets.push({
            type: 'line',
            label: `Bond ${index + 1}`,
            data: [
                { x: bond.from, y: 0.2 },  // Offset to show above the amino acids
                { x: bond.to, y: 0.2 }
            ],
            borderColor: 'blue',
            borderWidth: 2,
            pointRadius: 0,
            tension: 0.1
        });
    });
    
    // Create and render the chart
    if (window.structureChart) {
        window.structureChart.destroy();
    }
    
    window.structureChart = new Chart(ctx, {
        type: 'scatter',
        data: data,
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    display: false
                },
                x: {
                    ticks: {
                        callback: function(value) {
                            return sequence[value];
                        }
                    }
                }
            },
            plugins: {
                legend: {
                    display: false
                }
            }
        }
    });
}

// Initialize application when page loads
window.onload = main;