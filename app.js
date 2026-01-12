const state = {
  selectedCountries: new Set(),
  focusedCountry: null,
  countryFilter: null,
};

const targetCol = "Thriving Index";
const featureDefs = [
  { key: "help", label: "Help", col: "Count On to Help - Yes" },
  { key: "peace", label: "Peace", col: "Feel at Peace With Life - Yes" },
  { key: "balance", label: "Balance", col: "Various Aspects of Life in Balance - Yes" },
  { key: "health", label: "Health problems", col: "Health Problems - Yes" },
  { key: "calm", label: "Calm yesterday", col: "Experience Calmness Yesterday - Yes" },
  { key: "enjoy", label: "Enjoy work", col: "Enjoy the Work You Do Every Day - Yes" },
  { key: "choices", label: "Choices", col: "Many Choices in Type of Work - Yes" },
  { key: "exciting", label: "Exciting life", col: "Exciting Life or Calm Life - An exciting life" },
];

const purposeDefs = [
  {
    key: "purpose_good",
    label: "Be good at what you do",
    col: "Main Purpose in Life - Being good at what you do in your daily life",
  },
  {
    key: "purpose_family",
    label: "Caring for family",
    col: "Main Purpose in Life - Caring for family and close friends",
  },
  {
    key: "purpose_help",
    label: "Helping others",
    col: "Main Purpose in Life - Helping other people who need help",
  },
  { key: "purpose_none", label: "None of these", col: "Main Purpose in Life - (None of these)" },
  {
    key: "purpose_unknown",
    label: "Don't know",
    col: "Main Purpose in Life - (I do not know my purpose in life)",
  },
  { key: "purpose_dk", label: "DK", col: "Main Purpose in Life - (DK)" },
  { key: "purpose_refused", label: "Refused", col: "Main Purpose in Life - (Refused)" },
];

const purposeLegend = [
  { key: "purpose_good", label: "Be good at what you do", color: "#2c7fb8" },
  { key: "purpose_family", label: "Caring for family", color: "#7fcdbb" },
  { key: "purpose_help", label: "Helping others", color: "#41b6c4" },
  { key: "purpose_none", label: "None of these", color: "#f2c572" },
  { key: "purpose_unknown", label: "Don't know", color: "#e07a5f" },
  { key: "purpose_other", label: "DK/Refused", color: "#bcb8b1" },
];

const vibeDefs = [
  { key: "vibe_exciting", label: "Exciting", col: "Exciting Life or Calm Life - An exciting life" },
  { key: "vibe_calm", label: "Calm", col: "Exciting Life or Calm Life - A calm life" },
  { key: "vibe_both", label: "Both", col: "Exciting Life or Calm Life - (Both)" },
  { key: "vibe_neither", label: "Neither", col: "Exciting Life or Calm Life - (Neither)" },
  { key: "vibe_dk", label: "DK", col: "Exciting Life or Calm Life - (DK)" },
  { key: "vibe_refused", label: "Refused", col: "Exciting Life or Calm Life - (Refused)" },
];

const vibeLegend = [
  { key: "vibe_exciting", label: "Exciting", color: "#f58518", dash: "" },
  { key: "vibe_calm", label: "Calm", color: "#4c78a8", dash: "" },
  { key: "vibe_both", label: "Both", color: "#ffbf79", dash: "4 3" },
  { key: "vibe_neither", label: "Neither", color: "#9ecae9", dash: "4 3" },
];

const quintileMap = new Map([
  ["Poorest 20%", 1],
  ["Second 20%", 2],
  ["Middle 20%", 3],
  ["Fourth 20%", 4],
  ["Richest 20%", 5],
]);

const tooltip = d3.select("#tooltip");

function parsePercent(value) {
  if (value === undefined || value === null) {
    return NaN;
  }
  const cleaned = value.toString().trim();
  if (!cleaned) {
    return NaN;
  }
  return Number(cleaned.replace("%", "").replace(/\s/g, "").replace(",", "."));
}

function normalizeLabel(value) {
  return value.replace(/\u00a0/g, " ").replace(/\s+/g, " ").trim();
}

function buildColumns(header1, header2) {
  const filled = [];
  let current = "";
  header1.forEach((entry) => {
    if (entry.trim()) {
      current = entry.trim();
    }
    filled.push(current);
  });

  return header1.map((entry, i) => {
    if (i === 0) return "country";
    if (i === 1) return "group";
    if (i === 2) return "subgroup";
    const top = filled[i];
    const sub = header2[i].trim();
    if (sub && sub.toLowerCase() !== "rate") {
      return `${top} - ${sub}`;
    }
    return top;
  });
}

function transpose(matrix) {
  return matrix[0].map((_, i) => matrix.map((row) => row[i]));
}

function multiplyMatrices(a, b) {
  const result = Array.from({ length: a.length }, () => Array(b[0].length).fill(0));
  for (let i = 0; i < a.length; i += 1) {
    for (let k = 0; k < b.length; k += 1) {
      for (let j = 0; j < b[0].length; j += 1) {
        result[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  return result;
}

function multiplyMatrixVector(a, v) {
  return a.map((row) => row.reduce((sum, value, i) => sum + value * v[i], 0));
}

function invert(matrix) {
  const n = matrix.length;
  const augmented = matrix.map((row, i) => [
    ...row,
    ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)),
  ]);

  for (let i = 0; i < n; i += 1) {
    let pivot = i;
    for (let j = i + 1; j < n; j += 1) {
      if (Math.abs(augmented[j][i]) > Math.abs(augmented[pivot][i])) {
        pivot = j;
      }
    }
    if (pivot !== i) {
      [augmented[i], augmented[pivot]] = [augmented[pivot], augmented[i]];
    }
    const divisor = augmented[i][i];
    if (divisor === 0) {
      throw new Error("Matrix is singular.");
    }
    for (let j = 0; j < 2 * n; j += 1) {
      augmented[i][j] /= divisor;
    }
    for (let row = 0; row < n; row += 1) {
      if (row === i) continue;
      const factor = augmented[row][i];
      for (let col = 0; col < 2 * n; col += 1) {
        augmented[row][col] -= factor * augmented[i][col];
      }
    }
  }

  return augmented.map((row) => row.slice(n));
}

function linearRegression(X, y) {
  const X1 = X.map((row) => [1, ...row]);
  const Xt = transpose(X1);
  const XtX = multiplyMatrices(Xt, X1);
  const XtXInv = invert(XtX);
  const XtY = multiplyMatrixVector(Xt, y);
  const beta = multiplyMatrixVector(XtXInv, XtY);

  const predictions = X1.map((row) =>
    row.reduce((sum, value, i) => sum + value * beta[i], 0),
  );
  const yMean = d3.mean(y);
  const ssTot = d3.sum(y, (value) => (value - yMean) ** 2);
  const ssRes = d3.sum(predictions, (pred, i) => (y[i] - pred) ** 2);
  const r2 = ssTot ? 1 - ssRes / ssTot : 0;

  return {
    intercept: beta[0],
    coefficients: beta.slice(1),
    r2,
    predictions,
  };
}

function summarizeCoefficients(coefData, r2, label) {
  const sorted = [...coefData].sort((a, b) => b.coef - a.coef);
  const topPositive = sorted.slice(0, 2).map((d) => d.label).join(", ");
  const topNegative = sorted.slice(-2).map((d) => d.label).join(", ");
  const text =
    `${label}: R2 = ${r2.toFixed(2)}. ` +
    `Strongest positives: ${topPositive}. ` +
    `Strongest negatives: ${topNegative}.`;
  d3.select("#model-summary").text(text);
}

function fitModel(rows) {
  const X = rows.map((row) => featureDefs.map((def) => row[def.key]));
  const y = rows.map((row) => row.thriving);
  const means = featureDefs.map((_, i) => d3.mean(X, (row) => row[i]));
  const stds = featureDefs.map((_, i) => d3.deviation(X, (row) => row[i]) || 1);
  const Xscaled = X.map((row) => row.map((value, i) => (value - means[i]) / stds[i]));
  const regression = linearRegression(Xscaled, y);

  return {
    means,
    stds,
    intercept: regression.intercept,
    coefficients: regression.coefficients,
    r2: regression.r2,
  };
}

function predictRow(row, model) {
  const scaled = featureDefs.map(
    (def, i) => (row[def.key] - model.means[i]) / model.stds[i],
  );
  return model.intercept + d3.sum(scaled, (value, i) => value * model.coefficients[i]);
}

function incomeGradient(rows) {
  const data = rows
    .filter((d) => d.group === "Per Capita Income Quintiles" && d.quintile)
    .map((d) => ({ x: d.quintile, y: d.thriving }));
  if (data.length < 2) {
    return { slope: null, r2: null, gap: null };
  }

  const xs = data.map((d) => d.x);
  const ys = data.map((d) => d.y);
  const xMean = d3.mean(xs);
  const yMean = d3.mean(ys);
  const slope =
    d3.sum(xs, (x, i) => (x - xMean) * (ys[i] - yMean)) /
    d3.sum(xs, (x) => (x - xMean) ** 2);
  const intercept = yMean - slope * xMean;
  const ssTot = d3.sum(ys, (y) => (y - yMean) ** 2);
  const ssRes = d3.sum(ys, (y, i) => (y - (slope * xs[i] + intercept)) ** 2);
  const r2 = ssTot ? 1 - ssRes / ssTot : 0;

  const q1 = ys.filter((_, i) => xs[i] === 1);
  const q5 = ys.filter((_, i) => xs[i] === 5);
  const gap = q1.length && q5.length ? d3.mean(q5) - d3.mean(q1) : null;

  return { slope, r2, gap };
}

function loadData() {
  return d3.text("wrangle_dataset.csv").then((text) => {
    const rows = d3.dsvFormat(";").parseRows(text);
    const header1 = rows[0];
    const header2 = rows[1];
    const columns = buildColumns(header1, header2);
    const colIndex = new Map(columns.map((col, i) => [normalizeLabel(col), i]));
    const colIdx = (name) => colIndex.get(normalizeLabel(name));

    let currentCountry = "";
    let currentGroup = "";
    const records = [];

    rows.slice(3).forEach((row) => {
      const filled = row.slice();
      while (filled.length < columns.length) {
        filled.push("");
      }
      const country = (filled[0] || "").trim() || currentCountry;
      const group = (filled[1] || "").trim() || currentGroup;
      const subgroup = (filled[2] || "").trim();

      if ((filled[0] || "").trim()) currentCountry = country;
      if ((filled[1] || "").trim()) currentGroup = group;

      if (!country || !subgroup) return;

      const record = { country, group, subgroup };
      record.thriving = parsePercent(filled[colIdx(targetCol)]);

      featureDefs.forEach((def) => {
        record[def.key] = parsePercent(filled[colIdx(def.col)]);
      });

      purposeDefs.forEach((def) => {
        record[def.key] = parsePercent(filled[colIdx(def.col)]);
      });

      vibeDefs.forEach((def) => {
        record[def.key] = parsePercent(filled[colIdx(def.col)]);
      });

      record.quintile = quintileMap.get(subgroup) || null;
      if (Number.isFinite(record.thriving)) {
        record.rowId = records.length;
        records.push(record);
      }
    });

    return records;
  });
}

function buildCharts(records) {
  const regressionRows = records.filter((row) =>
    featureDefs.every((def) => Number.isFinite(row[def.key])) &&
      Number.isFinite(row.thriving),
  );

  const baseModel = fitModel(regressionRows);

  const byCountry = d3.group(regressionRows, (d) => d.country);
  const countryStats = Array.from(byCountry, ([country, rows]) => {
    const gradient = incomeGradient(rows);
    return {
      country,
      avgThriving: d3.mean(rows, (d) => d.thriving),
      avgBalance: d3.mean(rows, (d) => d.balance),
      avgPeace: d3.mean(rows, (d) => d.peace),
      incomeSlope: gradient.slope,
      incomeR2: gradient.r2,
      incomeGap: gradient.gap,
    };
  });

  setupScatter(countryStats);
  setupProfileChart();
  setupFitChart();
  setupPurposeChart();
  setupVibeChart();

  state.data = {
    countryStats,
    rows: records,
    regressionRows,
  };
  state.modelBase = baseModel;

  initCountryFilters(records);
  updateAll();

  d3.select("#clear-selection").on("click", () => {
    state.selectedCountries.clear();
    state.focusedCountry = null;
    if (state.scatter?.brushG) {
      state.scatter.brushG.call(state.scatter.brush.move, null);
    }
    updateAll();
  });
}

function setupScatter(countryStats) {
  const container = d3.select("#scatter-chart").node();
  const width = container.clientWidth;
  const height = container.clientHeight;
  const margin = { top: 20, right: 20, bottom: 50, left: 60 };

  const svg = d3.select("#scatter").attr("viewBox", `0 0 ${width} ${height}`);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;

  const data = countryStats.filter((d) => Number.isFinite(d.incomeSlope));
  const xScale = d3
    .scaleLinear()
    .domain(d3.extent(data, (d) => d.incomeSlope))
    .nice()
    .range([0, innerWidth]);
  const yScale = d3
    .scaleLinear()
    .domain(d3.extent(data, (d) => d.avgThriving))
    .nice()
    .range([innerHeight, 0]);
  const colorScale = d3
    .scaleSequential()
    .domain(d3.extent(data, (d) => d.avgBalance))
    .interpolator(d3.interpolateYlGnBu);
  const sizeScale = d3
    .scaleSqrt()
    .domain(d3.extent(data, (d) => d.incomeGap || 0))
    .range([4, 14]);

  g.append("g")
    .attr("class", "axis")
    .attr("transform", `translate(0,${innerHeight})`)
    .call(d3.axisBottom(xScale).ticks(6));

  g.append("g").attr("class", "axis").call(d3.axisLeft(yScale).ticks(6));

  g.append("text")
    .attr("x", innerWidth / 2)
    .attr("y", innerHeight + 40)
    .attr("text-anchor", "middle")
    .attr("fill", "#5c6a60")
    .text("Income gradient (thriving change per quintile)");

  g.append("text")
    .attr("x", -innerHeight / 2)
    .attr("y", -44)
    .attr("transform", "rotate(-90)")
    .attr("text-anchor", "middle")
    .attr("fill", "#5c6a60")
    .text("Average thriving (%)");

  const points = g
    .append("g")
    .attr("class", "points")
    .selectAll("circle")
    .data(data, (d) => d.country)
    .join("circle")
    .attr("class", "point")
    .attr("cx", (d) => xScale(d.incomeSlope))
    .attr("cy", (d) => yScale(d.avgThriving))
    .attr("r", (d) => sizeScale(d.incomeGap || 0))
    .attr("fill", (d) => colorScale(d.avgBalance))
    .on("mouseover", (event, d) => showTooltip(event, d))
    .on("mousemove", (event) => moveTooltip(event))
    .on("mouseout", hideTooltip)
    .on("click", (_, d) => {
      state.focusedCountry = d.country;
      updateAll();
    });

  const brush = d3
    .brush()
    .extent([
      [0, 0],
      [innerWidth, innerHeight],
    ])
    .on("brush end", (event) => {
      if (!event.selection) {
        state.selectedCountries.clear();
        updateAll();
        return;
      }
      const [[x0, y0], [x1, y1]] = event.selection;
      const selected = data.filter((d) => {
        const x = xScale(d.incomeSlope);
        const y = yScale(d.avgThriving);
        return x0 <= x && x <= x1 && y0 <= y && y <= y1;
      });
      state.selectedCountries = new Set(selected.map((d) => d.country));
      updateAll();
    });

  const brushG = g.append("g").attr("class", "scatter-brush").call(brush);

  const legend = d3.select("#scatter-legend");
  const colorDomain = colorScale.domain();
  const gapDomain = sizeScale.domain();
  legend.html(
    `<div>Balance: ${colorDomain[0].toFixed(1)}%</div>` +
      `<div class="bar"></div>` +
      `<div>${colorDomain[1].toFixed(1)}%</div>` +
      `<div>Gap size: ${gapDomain[0].toFixed(1)}-${gapDomain[1].toFixed(1)} pp</div>`,
  );

  state.scatter = { points, colorScale, sizeScale, brush, brushG };
}

function setupProfileChart() {
  const container = d3.select("#profile-chart").node();
  const width = container.clientWidth;
  const height = container.clientHeight;
  const margin = { top: 10, right: 20, bottom: 35, left: 140 };

  const svg = d3.select("#profile").attr("viewBox", `0 0 ${width} ${height}`);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;

  g.append("g").attr("class", "axis x-axis");
  g.append("g").attr("class", "axis y-axis");

  g.append("text")
    .attr("x", innerWidth / 2)
    .attr("y", innerHeight + 28)
    .attr("text-anchor", "middle")
    .attr("fill", "#5c6a60")
    .text("Standardized coefficient");

  state.profile = { svg, g, innerWidth, innerHeight };
}

function setupFitChart() {
  const container = d3.select("#fit-chart").node();
  const width = container.clientWidth;
  const height = container.clientHeight;
  const margin = { top: 10, right: 20, bottom: 35, left: 60 };

  const svg = d3.select("#fit").attr("viewBox", `0 0 ${width} ${height}`);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;

  g.append("g").attr("class", "axis x-axis");
  g.append("g").attr("class", "axis y-axis");

  g.append("text")
    .attr("x", innerWidth / 2)
    .attr("y", innerHeight + 28)
    .attr("text-anchor", "middle")
    .attr("fill", "#5c6a60")
    .text("Predicted thriving (%)");

  g.append("text")
    .attr("x", -innerHeight / 2)
    .attr("y", -44)
    .attr("transform", "rotate(-90)")
    .attr("text-anchor", "middle")
    .attr("fill", "#5c6a60")
    .text("Actual thriving (%)");

  state.fit = { svg, g, innerWidth, innerHeight };
}

function setupPurposeChart() {
  const container = d3.select("#purpose-chart").node();
  const width = container.clientWidth;
  const height = container.clientHeight;
  const margin = { top: 10, right: 20, bottom: 28, left: 140 };

  const svg = d3.select("#purpose").attr("viewBox", `0 0 ${width} ${height}`);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;

  g.append("g").attr("class", "axis x-axis");
  g.append("g").attr("class", "axis y-axis");

  g.append("text")
    .attr("x", innerWidth / 2)
    .attr("y", innerHeight + 24)
    .attr("text-anchor", "middle")
    .attr("fill", "#5c6a60")
    .text("Income quintile (1 = poorest, 5 = richest)");

  g.append("text")
    .attr("x", -innerHeight / 2)
    .attr("y", -44)
    .attr("transform", "rotate(-90)")
    .attr("text-anchor", "middle")
    .attr("fill", "#5c6a60")
    .text("Share of responses (%)");

  renderLegend("#purpose-legend", purposeLegend);
  state.purpose = { svg, g, innerWidth, innerHeight };
}

function setupVibeChart() {
  const container = d3.select("#vibe-chart").node();
  const width = container.clientWidth;
  const height = container.clientHeight;
  const margin = { top: 10, right: 20, bottom: 28, left: 110 };

  const svg = d3.select("#vibe").attr("viewBox", `0 0 ${width} ${height}`);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;

  g.append("g").attr("class", "axis x-axis");
  g.append("g").attr("class", "axis y-axis");

  g.append("text")
    .attr("x", innerWidth / 2)
    .attr("y", innerHeight + 24)
    .attr("text-anchor", "middle")
    .attr("fill", "#5c6a60")
    .text("Income quintile (1 = poorest, 5 = richest)");

  g.append("text")
    .attr("x", -innerHeight / 2)
    .attr("y", -44)
    .attr("transform", "rotate(-90)")
    .attr("text-anchor", "middle")
    .attr("fill", "#5c6a60")
    .text("Share of responses (%)");

  renderLegend("#vibe-legend", vibeLegend);
  state.vibe = { svg, g, innerWidth, innerHeight };
}

function renderLegend(containerId, items) {
  const container = d3.select(containerId);
  container.html("");
  const entry = container
    .selectAll("div.legend-item")
    .data(items)
    .enter()
    .append("div")
    .attr("class", "legend-item");

  entry
    .append("span")
    .attr("class", "legend-swatch")
    .style("background", (d) => d.color);

  entry.append("span").text((d) => d.label);
}

function getActiveSelection() {
  if (state.focusedCountry) {
    return { set: new Set([state.focusedCountry]), label: state.focusedCountry };
  }
  if (state.selectedCountries.size) {
    return { set: state.selectedCountries, label: "Selected countries" };
  }
  return null;
}

function getPurposeSelection() {
  if (state.countryFilter) {
    return { set: new Set([state.countryFilter]), label: state.countryFilter };
  }
  return getActiveSelection();
}

function initCountryFilters(records) {
  const countries = Array.from(new Set(records.map((d) => d.country))).sort();
  const options = ["All countries", ...countries];

  const selects = d3.selectAll(".country-filter");
  selects.selectAll("option").remove();
  selects
    .selectAll("option")
    .data(options)
    .enter()
    .append("option")
    .attr("value", (d) => d)
    .text((d) => d);

  const current = state.countryFilter || "All countries";
  selects.property("value", current);

  selects.on("change", (event) => {
    const value = event.target.value;
    state.countryFilter = value === "All countries" ? null : value;
    d3.selectAll(".country-filter").property("value", value);
    updateAll();
  });
}

function buildPurposeValues(rows) {
  const values = {};
  purposeDefs.forEach((def) => {
    values[def.key] =
      d3.mean(rows, (d) => (Number.isFinite(d[def.key]) ? d[def.key] : undefined)) || 0;
  });
  values.purpose_other = (values.purpose_dk || 0) + (values.purpose_refused || 0);
  return {
    purpose_good: values.purpose_good,
    purpose_family: values.purpose_family,
    purpose_help: values.purpose_help,
    purpose_none: values.purpose_none,
    purpose_unknown: values.purpose_unknown,
    purpose_other: values.purpose_other,
  };
}


function updateScatter(countryStats) {
  if (!state.scatter) return;
  const { points } = state.scatter;
  points
    .classed("selected", (d) => state.selectedCountries.has(d.country))
    .classed("focused", (d) => d.country === state.focusedCountry);
}

function updateProfile(coefData) {
  const chart = state.profile;
  if (!chart) return;

  const data = [...coefData].sort((a, b) => a.coef - b.coef);
  const maxAbs = d3.max(data, (d) => Math.abs(d.coef)) || 1;
  const xScale = d3
    .scaleLinear()
    .domain([-maxAbs, maxAbs])
    .range([0, chart.innerWidth])
    .nice();
  const yScale = d3
    .scaleBand()
    .domain(data.map((d) => d.label))
    .range([chart.innerHeight, 0])
    .padding(0.25);

  chart.g.select(".x-axis").attr("transform", `translate(0,${chart.innerHeight})`).call(
    d3.axisBottom(xScale).ticks(5),
  );
  chart.g.select(".y-axis").call(d3.axisLeft(yScale));

  const zeroLine = chart.g.selectAll("line.zero").data([0]);
  zeroLine
    .join("line")
    .attr("class", "zero")
    .attr("x1", xScale(0))
    .attr("x2", xScale(0))
    .attr("y1", 0)
    .attr("y2", chart.innerHeight)
    .attr("stroke", "#b7c3b8")
    .attr("stroke-dasharray", "3 3");

  const bars = chart.g.selectAll("rect.coef").data(data, (d) => d.label);
  bars
    .join("rect")
    .attr("class", (d) =>
      d.coef >= 0 ? "coef bar-positive" : "coef bar-negative",
    )
    .attr("x", (d) => xScale(Math.min(0, d.coef)))
    .attr("y", (d) => yScale(d.label))
    .attr("height", yScale.bandwidth())
    .attr("width", (d) => Math.abs(xScale(d.coef) - xScale(0)));
}

function updatePurpose(rows) {
  const chart = state.purpose;
  if (!chart) return;

  const active = getPurposeSelection();
  const incomeRows = rows.filter(
    (row) => row.group === "Per Capita Income Quintiles" && row.quintile,
  );
  const selectedRows = active
    ? incomeRows.filter((row) => active.set.has(row.country))
    : incomeRows;

  const data = d3.range(1, 6).map((q) => {
    const subset = selectedRows.filter((row) => row.quintile === q);
    return { quintile: q, ...buildPurposeValues(subset) };
  });

  const keys = purposeLegend.map((d) => d.key);
  const stacked = d3.stack().keys(keys)(data);

  const xScale = d3
    .scaleBand()
    .domain(data.map((d) => d.quintile))
    .range([0, chart.innerWidth])
    .padding(0.25);
  const yScale = d3.scaleLinear().domain([0, 100]).range([chart.innerHeight, 0]).nice();

  chart.g.select(".x-axis").attr("transform", `translate(0,${chart.innerHeight})`).call(
    d3.axisBottom(xScale).tickFormat(d3.format("d")),
  );
  chart.g.select(".y-axis").call(d3.axisLeft(yScale).ticks(5));

  const groups = chart.g.selectAll("g.stack").data(stacked, (d) => d.key);
  groups
    .join("g")
    .attr("class", "stack")
    .attr("fill", (d) => purposeLegend.find((p) => p.key === d.key).color)
    .selectAll("rect")
    .data((d) =>
      d.map((segment) => ({
        key: d.key,
        quintile: segment.data.quintile,
        y0: segment[0],
        y1: segment[1],
      })),
    )
    .join("rect")
    .attr("x", (d) => xScale(d.quintile))
    .attr("width", xScale.bandwidth())
    .attr("y", (d) => yScale(d.y1))
    .attr("height", (d) => yScale(d.y0) - yScale(d.y1));

  const labelText = chart.g.selectAll("text.series-label").data([active ? active.label : "All countries"]);
  labelText
    .join("text")
    .attr("class", "series-label")
    .attr("x", chart.innerWidth - 6)
    .attr("y", 12)
    .attr("text-anchor", "end")
    .attr("fill", "#5c6a60")
    .text(active ? active.label : "All countries");
}

function updateVibe(rows) {
  const chart = state.vibe;
  if (!chart) return;

  const active = getPurposeSelection();
  const incomeRows = rows.filter(
    (row) => row.group === "Per Capita Income Quintiles" && row.quintile,
  );
  const selectedRows = active
    ? incomeRows.filter((row) => active.set.has(row.country))
    : incomeRows;

  const data = d3.range(1, 6).map((q) => {
    const subset = selectedRows.filter((row) => row.quintile === q);
    return {
      quintile: q,
      exciting:
        d3.mean(subset, (d) => (Number.isFinite(d.vibe_exciting) ? d.vibe_exciting : undefined)) ||
        0,
      calm:
        d3.mean(subset, (d) => (Number.isFinite(d.vibe_calm) ? d.vibe_calm : undefined)) || 0,
      both:
        d3.mean(subset, (d) => (Number.isFinite(d.vibe_both) ? d.vibe_both : undefined)) || 0,
      neither:
        d3.mean(subset, (d) => (Number.isFinite(d.vibe_neither) ? d.vibe_neither : undefined)) ||
        0,
    };
  });

  const xScale = d3.scaleLinear().domain([1, 5]).range([0, chart.innerWidth]);
  const yScale = d3.scaleLinear().domain([0, 100]).range([chart.innerHeight, 0]).nice();

  chart.g.select(".x-axis").attr("transform", `translate(0,${chart.innerHeight})`).call(
    d3.axisBottom(xScale).ticks(5).tickFormat(d3.format("d")),
  );
  chart.g.select(".y-axis").call(d3.axisLeft(yScale).ticks(5));

  const series = [
    { key: "vibe_exciting", label: "Exciting", accessor: (d) => d.exciting },
    { key: "vibe_calm", label: "Calm", accessor: (d) => d.calm },
    { key: "vibe_both", label: "Both", accessor: (d) => d.both },
    { key: "vibe_neither", label: "Neither", accessor: (d) => d.neither },
  ];

  const line = d3
    .line()
    .defined((d) => Number.isFinite(d.value))
    .x((d) => xScale(d.quintile))
    .y((d) => yScale(d.value));

  const lines = chart.g.selectAll("path.vibe-line").data(series, (d) => d.key);
  lines
    .join("path")
    .attr("class", "vibe-line")
    .attr("fill", "none")
    .attr("stroke", (d) => vibeLegend.find((v) => v.key === d.key).color)
    .attr("stroke-width", 2)
    .attr("stroke-dasharray", (d) => vibeLegend.find((v) => v.key === d.key).dash)
    .attr("d", (d) =>
      line(data.map((row) => ({ quintile: row.quintile, value: d.accessor(row) }))),
    );

  const labelText = chart.g.selectAll("text.series-label").data([active ? active.label : "All countries"]);
  labelText
    .join("text")
    .attr("class", "series-label")
    .attr("x", chart.innerWidth - 6)
    .attr("y", 12)
    .attr("text-anchor", "end")
    .attr("fill", "#5c6a60")
    .text(active ? active.label : "All countries");
}

function updateFit(rows, model) {
  const chart = state.fit;
  if (!chart) return;

  let label = "All countries";
  let focusSet = null;
  let selectionSet = null;
  if (state.focusedCountry) {
    focusSet = new Set([state.focusedCountry]);
    label = state.focusedCountry;
  } else if (state.selectedCountries.size) {
    selectionSet = state.selectedCountries;
    label = "Selected countries";
  }

  const maxVal =
    d3.max(rows, (d) => Math.max(d.thriving, predictRow(d, model))) || 100;
  const xScale = d3.scaleLinear().domain([0, maxVal]).range([0, chart.innerWidth]).nice();
  const yScale = d3.scaleLinear().domain([0, maxVal]).range([chart.innerHeight, 0]).nice();

  chart.g.select(".x-axis").attr("transform", `translate(0,${chart.innerHeight})`).call(
    d3.axisBottom(xScale).ticks(5),
  );
  chart.g.select(".y-axis").call(d3.axisLeft(yScale).ticks(5));

  const ref = chart.g.selectAll("line.fit-line").data([0]);
  ref
    .join("line")
    .attr("class", "fit-line")
    .attr("x1", xScale(0))
    .attr("y1", yScale(0))
    .attr("x2", xScale(maxVal))
    .attr("y2", yScale(maxVal));

  const points = chart.g.selectAll("circle.fit-point").data(rows, (d) => d.rowId);
  points
    .join(
      (enter) => enter.append("circle").attr("class", "fit-point").attr("r", 3),
      (update) => update,
      (exit) => exit.remove(),
    )
    .attr("cx", (d) => xScale(predictRow(d, model)))
    .attr("cy", (d) => yScale(d.thriving))
    .classed(
      "dim",
      (d) => (focusSet || selectionSet) && !(focusSet || selectionSet).has(d.country),
    )
    .classed("focused", (d) => focusSet && focusSet.has(d.country))
    .classed("selected", (d) => selectionSet && selectionSet.has(d.country));

  const labelText = chart.g.selectAll("text.series-label").data([label]);
  labelText
    .join("text")
    .attr("class", "series-label")
    .attr("x", chart.innerWidth - 6)
    .attr("y", 12)
    .attr("text-anchor", "end")
    .attr("fill", "#5c6a60")
    .text(label);
}

function updateSelectionNote() {
  const count = state.selectedCountries.size;
  const label = count ? `Selected: ${count} countries` : "Selected: All countries";
  d3.select("#selection-note").text(label);
}

function updateAll() {
  const { countryStats, regressionRows, rows } = state.data;
  const active = getActiveSelection();
  const modelRows = active
    ? regressionRows.filter((row) => active.set.has(row.country))
    : regressionRows;

  let model = state.modelBase;
  let modelLabel = "All countries";
  if (active) {
    modelLabel = active.label;
    if (modelRows.length >= 20) {
      try {
        model = fitModel(modelRows);
      } catch (err) {
        model = state.modelBase;
        modelLabel = "All countries";
      }
    } else {
      modelLabel = "All countries";
    }
  }

  const coefData = featureDefs.map((def, i) => ({
    ...def,
    coef: model.coefficients[i],
  }));

  summarizeCoefficients(coefData, model.r2, modelLabel);
  updateSelectionNote();
  updateScatter(countryStats);
  updateProfile(coefData);
  updateFit(regressionRows, model);
  updatePurpose(rows);
  updateVibe(rows);
}

function showTooltip(event, d) {
  tooltip
    .classed("show", true)
    .style("left", `${event.pageX + 12}px`)
    .style("top", `${event.pageY - 24}px`)
    .html(
      `<strong>${d.country}</strong><br/>` +
        `Avg thriving: ${d.avgThriving.toFixed(1)}%<br/>` +
        `Income slope: ${d.incomeSlope.toFixed(2)}<br/>` +
        `Income gap: ${d.incomeGap ? d.incomeGap.toFixed(1) : "n/a"} pp`,
    );
}

function moveTooltip(event) {
  tooltip.style("left", `${event.pageX + 12}px`).style("top", `${event.pageY - 24}px`);
}

function hideTooltip() {
  tooltip.classed("show", false);
}

loadData().then(buildCharts).catch((error) => {
  // eslint-disable-next-line no-console
  console.error("Failed to load data:", error);
});
