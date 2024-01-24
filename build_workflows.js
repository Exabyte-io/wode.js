const fs = require("fs");
const path = require("path");
const yaml = require("js-yaml");

const allApplications = [
    "espresso",
    "jupyterLab",
    "exabyteml",
    "nwchem",
    "python",
    "shell",
    "vasp",
    "deepmd",
];

const allWorkflows = { workflows: {}, subworkflows: {} };

const JSONstringifyOrder = (obj, space) => {
    const allKeys = new Set();
    // eslint-disable-next-line no-sequences
    JSON.stringify(obj, (key, value) => (allKeys.add(key), value));
    return JSON.stringify(obj, Array.from(allKeys).sort(), space);
};

const loadFile = (name, dir, file, type) => {
    const obj = fs.readFileSync(path.resolve(dir, file), "utf8");
    const key = file.split(".")[0];
    allWorkflows[type][name][key] = yaml.load(obj);
};

allApplications.forEach((name) => {
    allWorkflows.workflows[name] = {};
    allWorkflows.subworkflows[name] = {};
    const wfDir = path.resolve(__dirname, "assets", "workflows", name);
    const swDir = path.resolve(__dirname, "assets", "subworkflows", name);
    try {
        const wfFiles = fs.readdirSync(wfDir);
        const swFiles = fs.readdirSync(swDir);
        console.log(
            `Building ${name}: ${wfFiles.length} workflow(s) and ${swFiles.length} subworkflow(s)`,
        );
        wfFiles.forEach((file) => loadFile(name, wfDir, file, "workflows"));
        swFiles.forEach((file) => loadFile(name, swDir, file, "subworkflows"));
    } catch (e) {
        console.log(e);
    }
});

const write_path = "workflows.js";
// write to src for unit test coverage simplicity
fs.writeFileSync(
    `./src/workflows/${write_path}`,
    "module.exports = {workflowData: " + JSONstringifyOrder(allWorkflows) + "}",
    "utf8",
);
fs.writeFileSync(
    `./dist/workflows/${write_path}`,
    "module.exports = {workflowData: " + JSONstringifyOrder(allWorkflows) + "}",
    "utf8",
);
