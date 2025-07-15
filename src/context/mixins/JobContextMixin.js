const defaultJob = {
    workflow: {
        subworkflows: [],
        units: [],
    },
    status: "pre-submission",
    compute: {
        queue: "D",
        nodes: 1,
        ppn: 1,
        timeLimit: "3600",
    },
    _project: {
        _id: "",
    },
};

export function jobContextMixin(item) {
    const properties = {
        isEdited: false,

        _job: defaultJob,

        get job() {
            return this._job;
        },

        initJobContextMixin() {
            const { config } = this;
            this._job = (config.context && config.context.job) || defaultJob;
            this.isEdited = false; // we always get the `defaultData` (recalculated from scratch, not persistent)
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
