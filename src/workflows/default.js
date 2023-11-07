export default {
    name: "New Workflow",
    properties: [],
    subworkflows: [
        {
            _id: "c6e9dbbee8929de01f4e76ee",
            application: {
                name: "espresso",
                summary: "Quantum Espresso",
                version: "6.3",
            },
            model: {
                method: {
                    subtype: "us",
                    type: "pseudopotential",
                },
                subtype: "gga",
                type: "dft",
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
