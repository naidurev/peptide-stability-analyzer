// js/main.js

// Global variable for Pyodide
let pyodide = null;

async function main() {
  try {
    // Show loading indicator and hide app container
    document.getElementById("loading-indicator").style.display = "flex";
    document.getElementById("app-container").style.display = "none";

    // Load Pyodide
    pyodide = await loadPyodide();
    console.log("Pyodide loaded successfully");

    // Load required packages (numpy and json are needed)
    await pyodide.loadPackagesFromImports(`
      import numpy
      import json
    `);

    // Fetch and run the Python code from the file
    const pythonCode = await fetchPythonCode();
    await pyodide.runPythonAsync(pythonCode);
    console.log("Python environment ready");

    // Hide loading indicator and show app container
    document.getElementById("loading-indicator").style.display = "none";
    document.getElementById("app-container").style.display = "block";

    // Set up UI event listeners
    setupEventListeners();
  } catch (error) {
    console.error("Error initializing application:", error);
    const loadingIndicator = document.getElementById("loading-indicator");
    if (loadingIndicator) {
      loadingIndicator.innerHTML = `
        <div class="alert alert-danger">
          <h4>Error Loading Application</h4>
          <p>${error.message}</p>
          <p>Please try refreshing the page or using a different browser.</p>
        </div>
      `;
    }
  }
}

async function fetchPythonCode() {
  try {
    const response = await fetch("python/peptide_analyzer.py");
    if (!response.ok) {
      throw new Error(`Failed to fetch Python code: ${response.status}`);
    }
    return await response.text();
  } catch (error) {
    console.error("Error fetching Python code:", error);
    throw error;
  }
}

function setupEventListeners() {
  const analyzeButton = document.getElementById("analyze-button");
  const sequenceInput = document.getElementById("sequence-input");

  analyzeButton.addEventListener("click", () => {
    analyzeSequence(sequenceInput.value);
  });

  sequenceInput.addEventListener("keypress", (event) => {
    if (event.key === "Enter") {
      analyzeSequence(sequenceInput.value);
    }
  });
}

async function analyzeSequence(sequence) {
  sequence = sequence.trim().toUpperCase();
  if (!sequence) {
    showStatusMessage("Please enter a peptide sequence.", "danger");
    return;
  }
  const validAminoAcids = "ACDEFGHIKLMNPQRSTVWY";
  if (![...sequence].every((aa) => validAminoAcids.includes(aa))) {
    showStatusMessage(
      "Invalid amino acid sequence. Please use standard one-letter codes.",
      "danger"
    );
    return;
  }

  try {
    // Run Python analysis code; ensure the last line returns "result"
    const result = await pyodide.runPythonAsync(`
import json
try:
    analyzer = PeptideStabilityAnalyzer("${sequence}")
    patterns = analyzer.get_possible_disulfide_patterns()
    if not patterns:
        result = json.dumps({"error": "No valid disulfide patterns found (need even number of cysteines)."})
    else:
        results = []
        for pattern in patterns:
            metrics = analyzer.evaluate_pattern_stability(pattern)
            results.append({
                "pattern": pattern,
                "metrics": metrics
            })
        results.sort(key=lambda x: x["metrics"]["total_energy"])
        result = json.dumps(results)
except Exception as e:
    result = json.dumps({"error": str(e)})
result
    `);

    const data = JSON.parse(result);
    if (data.error) {
      showStatusMessage(data.error, "danger");
      return;
    }
    displayResults(data, sequence);
    showStatusMessage("Analysis complete!", "success");
  } catch (error) {
    console.error("Error analyzing sequence:", error);
    showStatusMessage(`Error analyzing sequence: ${error.message}`, "danger");
  }
}

function showStatusMessage(message, type = "info") {
  const statusMessage = document.getElementById("status-message");
  if (!statusMessage) return;
  statusMessage.textContent = message;
  statusMessage.className = `alert alert-${type}`;
  statusMessage.style.display = "block";
  if (type === "success") {
    setTimeout(() => {
      statusMessage.style.display = "none";
    }, 5000);
  }
}

function displayResults(results, sequence) {
  // Display Energy Components in a bar chart
  createEnergyChart(results.slice(0, 5));
  // Populate the Pattern Analysis table
  createPatternTable(results);
  // Render sequence visualization and best pattern
  createSequenceVisualization(sequence, results[0].pattern);
}

function createEnergyChart(topResults) {
  const ctx = document.getElementById("energy-chart").getContext("2d");
  const labels = topResults.map((_, index) => `Pattern ${index + 1}`);
  const datasets = [];
  const componentKeys = Object.keys(topResults[0].metrics.energy_components);
  const colors = [
    "rgba(54, 162, 235, 0.7)",
    "rgba(255, 99, 132, 0.7)",
    "rgba(75, 192, 192, 0.7)",
    "rgba(255, 206, 86, 0.7)",
    "rgba(153, 102, 255, 0.7)",
    "rgba(255, 159, 64, 0.7)"
  ];
  componentKeys.forEach((component, index) => {
    const data = topResults.map(
      (result) => result.metrics.energy_components[component]
    );
    datasets.push({
      label: component.replace(/_/g, " "),
      data: data,
      backgroundColor: colors[index % colors.length],
      borderColor: colors[index % colors.length].replace("0.7", "1"),
      borderWidth: 1
    });
  });
  if (window.energyChart) {
    window.energyChart.destroy();
  }
  window.energyChart = new Chart(ctx, {
    type: "bar",
    data: {
      labels: labels,
      datasets: datasets
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        x: { stacked: true },
        y: {
          stacked: true,
          title: { display: true, text: "Energy (arbitrary units)" }
        }
      }
    }
  });
}

function createPatternTable(results) {
  const tableBody = document.getElementById("pattern-table-body");
  tableBody.innerHTML = "";
  results.forEach((result, index) => {
    const row = document.createElement("tr");
    // Pattern column
    const patternCell = document.createElement("td");
    const patternStr = result.pattern
      .map(([c1, c2]) => `C${c1 + 1}-C${c2 + 1}`)
      .join(" & ");
    patternCell.textContent = patternStr;
    row.appendChild(patternCell);
    // Total energy column
    const energyCell = document.createElement("td");
    energyCell.textContent = result.metrics.total_energy.toFixed(2);
    row.appendChild(energyCell);
    // Loop sizes column
    const loopCell = document.createElement("td");
    loopCell.textContent = result.metrics.loop_sizes.join(", ");
    row.appendChild(loopCell);
    // Highlight best pattern (first row)
    if (index === 0) {
      row.classList.add("table-success");
    }
    tableBody.appendChild(row);
  });
}

function createSequenceVisualization(sequence, bestPattern) {
  const sequenceDisplay = document.getElementById("sequence-display");
  sequenceDisplay.innerHTML = "";
  [...sequence].forEach((aa, index) => {
    const span = document.createElement("span");
    span.textContent = aa;
    span.classList.add("amino-acid");
    span.dataset.index = index;
    if (aa === "C") {
      span.classList.add("cysteine");
    }
    sequenceDisplay.appendChild(span);
  });
  createStructureChart(sequence, bestPattern);
}

function createStructureChart(sequence, bestPattern) {
  const ctx = document.getElementById("structure-chart").getContext("2d");
  const bondPairs = bestPattern.map(([c1, c2]) => ({ from: c1, to: c2 }));
  const aaPositions = [...sequence].map((_, index) => ({ x: index, y: 0 }));
  const data = {
    labels: [...sequence].map((aa, idx) => `${aa}${idx + 1}`),
    datasets: [
      {
        type: "scatter",
        label: "Amino Acids",
        data: aaPositions,
        backgroundColor: [...sequence].map((aa) =>
          aa === "C" ? "red" : "black"
        ),
        pointRadius: 8,
        pointHoverRadius: 10
      }
    ]
  };
  bondPairs.forEach((bond, index) => {
    data.datasets.push({
      type: "line",
      label: `Bond ${index + 1}`,
      data: [
        { x: bond.from, y: 0.2 },
        { x: bond.to, y: 0.2 }
      ],
      borderColor: "blue",
      borderWidth: 2,
      pointRadius: 0,
      tension: 0.1
    });
  });
  if (window.structureChart) {
    window.structureChart.destroy();
  }
  window.structureChart = new Chart(ctx, {
    type: "scatter",
    data: data,
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        y: { display: false },
        x: {
          ticks: {
            callback: function (value) {
              return sequence[value];
            }
          }
        }
      },
      plugins: {
        legend: { display: false }
      }
    }
  });
}

window.onload = main;

