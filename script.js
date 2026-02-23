const inputs = [
  document.getElementById("meanInput"),
  document.getElementById("standardDeviationInput"),
  document.getElementById("transformationFunctionInput"),
  document.getElementById("alphaInput"),
  document.getElementById("betaInput"),
  document.getElementById("kappaInput"),
];

let Traces = [];
//let mu = 1;
let standard_deviation = 1;
let transformation_function = "x";

let x_points;
let y_points;

let old_mean = 0;
let old_std_dev = 0;
let new_mean = 0;
let new_std_dev = 0;
let new_covariance = 0;

let n;
let alpha;
let beta;
let kappa;
let lambda;

let sigmas = [];

let weights_mean = [];
let weights_covariance = [];

let transformedPoints = [];
//let mu = 1;
//let standard_deviation = 1;

//    let gaussian_expression = "(1 / (standard_deviation * Math.sqrt(2 * Math.PI))) * Math.exp(-Math.pow(x - mu, 2) / (2 * standard_deviation * standard_deviation));"

function generatePoints(expression, x_min, x_max, step) {
  const x_points = [];
  const y_points = [];

  for (let x = x_min; x <= x_max; x += step) {
    x_points.push(x);
    const val = eval(expression);
    y_points.push(val);
  }

  return [x_points, y_points]; // vrací pole
}
function drawCurve(x_points, y_points, name, color) {
  Traces.push({
    x: x_points,
    y: y_points,
    mode: "lines",
    line: { color: color },
    name: name,
  });
}
// Standard Normal variate using Box-Muller transform.
function gaussianRandom(mean, stdev) {
  const u = 1 - Math.random(); // Converting [0,1) to (0,1]
  const v = Math.random();
  const z = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
  // Transform to the desired mean and standard deviation:
  return z * stdev + mean;
}
function transformPoints(expression, x_points) {
  for (let i = 0; i < x_points.length; i++) {
    let x = x_points[i];
    x_points[i] = eval(expression);
  }
  return x_points;
}
function generateGaussian(mean, std_dev, min_x, max_x, step) {
  const x_points = [];
  const y_points = [];

  for (let x = min_x; x <= max_x; x += step) {
    x_points.push(x);
    const y =
      (1 / (std_dev * Math.sqrt(2 * Math.PI))) *
      Math.exp(-Math.pow(x - mean, 2) / (2 * std_dev * std_dev));
    y_points.push(y);
  }

  return [x_points, y_points];
}
function graphWeightedPoints(
  x_points,
  y_points,
  weights,
  name,
  colorPositive,
  colorNegative
) {
  const totalAbs = weights.reduce((sum, w) => sum + Math.abs(w), 0);

  Traces.push({
    x: x_points,
    y: y_points,
    mode: "markers",
    type: "scatter",
    marker: {
      color: weights.map((w) => (w >= 0 ? colorPositive : colorNegative)),
      size: weights.map((w) => (Math.abs(w) / totalAbs) * 30),
      symbol: "circle",
    },
    name: name,
  });
}
function graphPoints(transformedPoints, y_points, name, color) {
  Traces.push({
    x: transformedPoints,
    y: y_points,
    mode: "markers",
    type: "scatter",
    marker: {
      color: color,
      size: 10,
      symbol: "circle",
    },
    name: name,
  });
}

function unscentedTransform(mean, standard_deviation, transformation_function) {
  sigmas[0] = mean;
  // generate sigma points
  // only 1D; L = 2n + 1
  n = 1;
  alpha = parseFloat(document.getElementById("alphaInput").value);
  beta = parseFloat(document.getElementById("betaInput").value);
  kappa = parseFloat(document.getElementById("kappaInput").value);

  lambda = alpha ** 2 * (n + kappa) - n;

  // Numerical stability check
  if (Math.abs(n + lambda) < 1e-10) {
    console.error(
      "Numerically instable parameters! Try larger alpha or kappa."
    );
    return [mean, standard_deviation ** 2];
  }

  let spread = Math.sqrt((n + lambda) * standard_deviation ** 2);
  sigmas[1] = mean + spread;
  sigmas[2] = mean - spread;

  // generate weights
  weights_mean[0] = lambda / (n + lambda);
  weights_covariance[0] = lambda / (n + lambda) + (1 - alpha ** 2 + beta);

  for (let i = 1; i <= 2 * n; i++) {
    let tmp = 1 / (2 * (n + lambda));
    weights_mean[i] = tmp;
    weights_covariance[i] = tmp;
  }

  //  console.log("alpha:", alpha, "beta:", beta, "kappa:", kappa);
  //  console.log("lambda:", lambda);
  //  console.log("n + lambda:", n + lambda);
  //  console.log("weights_mean:", weights_mean);
  //  console.log("weights_covariance:", weights_covariance[0]);

  graphWeightedPoints(
    sigmas,
    [-0.1, -0.1, -0.1],
    weights_mean,
    "Sigma Points",
    "blue",
    "red"
  );

  let transformedSigmaPoints = [];
  let x = 0;
  for (let i = 0; i < 3; i++) {
    x = sigmas[i];
    transformedSigmaPoints[i] = eval(transformation_function);
  }

  graphWeightedPoints(
    transformedSigmaPoints,
    [-0.3, -0.3, -0.3],
    weights_mean,
    "Transformed Sigma Points",
    "blue",
    "red"
  );

  let new_mean = 0;
  let new_covariance = 0;
  for (let i = 0; i < 2 * n + 1; i++) {
    new_mean += weights_mean[i] * transformedSigmaPoints[i];
  }
  for (let i = 0; i < 2 * n + 1; i++) {
    new_covariance +=
      weights_covariance[i] * (transformedSigmaPoints[i] - new_mean) ** 2;
  }

  return [new_mean, new_covariance];
}

function updateLatex() {
  const solutionString = `\\mu = ${old_mean.toFixed(
    3
  )}, \\quad \\sigma = ${old_std_dev.toFixed(
    3
  )}, \\quad \\mu' = ${new_mean.toFixed(
    3
  )}, \\quad \\sigma' = ${new_std_dev.toFixed(3)}`;
  katex.render(solutionString, document.getElementById("solutionLatex"), {
    throwOnError: false,
  });

  const computationString = `\\lambda = ${alpha}^2 \\times (${n} + ${kappa}) - ${n} = ${lambda.toFixed(
    3
  )}`;
  katex.render(computationString, document.getElementById("lambdaLatex"), {
    throwOnError: false,
  });

  const sigmasString = `\\chi^0 = ${sigmas[0].toFixed(
    3
  )}, \\chi^1 = ${sigmas[1].toFixed(3)}, \\chi^2 = ${sigmas[2].toFixed(3)}`;
  katex.render(sigmasString, document.getElementById("sigmasLatex"), {
    throwOnError: false,
  });

  const weightsMeanString = `w^0_m = ${weights_mean[0].toFixed(
    3
  )}, w^1_m = ${weights_mean[1].toFixed(3)}, w^2_m = ${weights_mean[2].toFixed(
    3
  )}`;
  katex.render(weightsMeanString, document.getElementById("weightsMeanLatex"), {
    throwOnError: false,
  });

  const weightsCovarianceString = `w^0_c = ${weights_covariance[0].toFixed(
    3
  )}, w^1_c = ${weights_covariance[1].toFixed(
    3
  )}, w^2_c = ${weights_covariance[2].toFixed(3)}`;
  katex.render(
    weightsCovarianceString,
    document.getElementById("weightsCovarianceLatex"),
    {
      throwOnError: false,
    }
  );
  //newMeanComputationLatex
  const newMeanComputationString = `\\mu' = \\sum_{i=0}^{2n} w_m^{[i]} \\mathbf{y}^{[i]} = ${weights_mean[0].toFixed(
    3
  )} \\times ${transformedPoints[0].toFixed(3)} + ${weights_mean[1].toFixed(
    3
  )} \\times ${transformedPoints[1].toFixed(3)} + ${weights_mean[2].toFixed(
    3
  )} \\times ${transformedPoints[2].toFixed(3)} = ${new_mean.toFixed(3)}`;
  katex.render(
    newMeanComputationString,
    document.getElementById("newMeanComputationLatex"),
    {
      throwOnError: false,
    }
  );

  const newCovarianceComputationString = `\\sigma' = \\sum_{i=0}^{2n} w_c^{[i]} \\mathbf{y}^{[i]} = ${weights_covariance[0].toFixed(
    3
  )} \\times ${transformedPoints[0].toFixed(
    3
  )} + ${weights_covariance[1].toFixed(
    3
  )} \\times ${transformedPoints[1].toFixed(
    3
  )} + ${weights_covariance[2].toFixed(
    3
  )} \\times ${transformedPoints[2].toFixed(3)} = ${new_covariance.toFixed(3)}`;
  katex.render(
    newCovarianceComputationString,
    document.getElementById("newCovarianceComputationLatex"),
    {
      throwOnError: false,
    }
  );
}

function updatePlot() {
  let visibilityState = {};
  const plotDiv = document.getElementById("plot");
  if (plotDiv && plotDiv.data) {
    plotDiv.data.forEach(trace => {
      if (trace.name) {
        visibilityState[trace.name] = trace.visible; 
      }
    });
  }

  Traces = [];
  old_mean = parseFloat(document.getElementById("meanInput").value);
  old_std_dev = parseFloat(document.getElementById("standardDeviationInput").value);
  transformation_function = document.getElementById("transformationFunctionInput").value;

  const graphLayout = {
    title: "Unscented Transform",
    xaxis: { title: "x", range: [-2, 2] },
    yaxis: {
      title: "y",
      range: [-1.5, 1.5],
    },
  };

  [x_points, y_points] = generateGaussian(old_mean, old_std_dev, -2, 2, 0.01);
  drawCurve(x_points, y_points, "Original Gaussian", "green");

  [x_points, y_points] = generatePoints(transformation_function, -2, 2, 0.01);
  drawCurve(x_points, y_points, "Transformation Function", "red");

  // histogram
  var histogram_x = [];
  for (var i = 0; i < 50000; i++) {
    histogram_x[i] = gaussianRandom(old_mean, old_std_dev); 
  }
  histogram_x = transformPoints(transformation_function, histogram_x);
  Traces.push({
    x: histogram_x,
    type: "histogram",
    histnorm: "probability density",
    marker: { color: "rgba(0,128,128,0.5)" },
    name: "Monte Carlo Approximation",
    opacity: 0.5,
  });

  [new_mean, new_covariance] = unscentedTransform(
    old_mean,
    old_std_dev,
    transformation_function
  );

  new_std_dev = Math.sqrt(new_covariance);
  transformedPoints = [
    new_mean,
    new_mean - new_std_dev,
    new_mean + new_std_dev,
  ];

  graphPoints(
    transformedPoints,
    [-0.5, -0.5, -0.5],
    "Unscented Transformed",
    "orange"
  );

  [x_points, y_points] = generateGaussian(new_mean, new_std_dev, -5, 5, 0.01);
  drawCurve(x_points, y_points, "Transformed Gaussian", "blue");

  Traces.forEach(trace => {
    if (trace.name && visibilityState[trace.name] !== undefined) {
      trace.visible = visibilityState[trace.name];
    }
  });

  //Plotly.react("plot", Traces, graphLayout, {displayModeBar: false});
  Plotly.react("plot", Traces, graphLayout, {displayModeBar: false, responsive: true});
  updateLatex();
}

updatePlot();

inputs.forEach((input) => {
  input.addEventListener("input", updatePlot);
});
