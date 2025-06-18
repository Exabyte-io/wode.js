// import type { ContextProvider } from "@mat3ra/code/dist/js/context";
// import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
// import type { JobSchema } from "@mat3ra/esse/dist/js/types";

// type Job = JobSchema & {
//     parentJob?: Job;
// };

// type JobConfig = {
//     context?: {
//         job?: Job;
//     };
// };

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

// export type JobContextMixinType = {
//     isEdited?: boolean;
//     job: Job;
//     _job: Job;
//     initJobContextMixin: () => void;
// };

// export function jobContextMixin(item: ContextProvider) {
//     const properties = {
//         isEdited: false,

//         _job: defaultJob,

//         get job() {
//             return this._job;
//         },

//         initJobContextMixin() {
//             const config = this.config as JobConfig;
//             this._job = (config.context && config.context.job) || defaultJob;
//             this.isEdited = false; // we always get the `defaultData` (recalculated from scratch, not persistent)
//         },
//     } as JobContextMixinType & typeof item;

//     Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
// }

// export function JobContextMixin<T extends Constructor<ContextProvider>>(superclass: T) {
//     jobContextMixin(superclass.prototype);
//     return superclass as T & Constructor<JobContextMixinType>;
// }

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
