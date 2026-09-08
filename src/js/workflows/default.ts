import type { WorkflowSchema } from "./types";

const defaultWorkflowConfig: WorkflowSchema = {
    name: "New Workflow",
    properties: [],
    subworkflows: [
        {
            _id: "c6e9dbbee8929de01f4e76ee",
            application: {
                name: "espresso",
                shortName: "qe",
                summary: "Quantum ESPRESSO",
                build: "GNU",
                version: "6.3",
            },
            model: {
                method: {
                    subtype: "us",
                    type: "pseudopotential",
                },
                subtype: "gga",
                type: "dft",
                functional: "other",
            },
            name: "New Subworkflow",
            properties: [],
            units: [],
        },
    ],
    workflows: [],
    units: [
        {
            _id: "c6e9dbbee8929de01f4e76ee",
            flowchartId: "da2c090ede4dc2fa6e66647f",
            head: true,
            monitors: [],
            postProcessors: [],
            preProcessors: [],
            results: [],
            type: "subworkflow",
            name: "New Subworkflow",
        },
    ],
};

export default defaultWorkflowConfig;
