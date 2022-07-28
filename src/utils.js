import lodash from "lodash";

/**
 * @summary set the head of an array of units
 * @param units
 * @returns {Unit[]}
 */
export function setUnitsHead(units) {
    if (units.length > 0) {
        units[0].head = true;
        lodash.tail(units).map(x => x.head = false);
    }
    return units;

}

// TODO: fix setNextLinks on unit removal and convergence logic.
/**
 * @summary Re-establishes the linked `next => flowchartId` logic in an array of units
 * @params units {Unit[]}
 * @returns units {Unit[]}
 */
export function setNextLinks(units) {
    const flowchartIds = units.map(u => u.flowchartId);
    for (let i = 0; i < units.length - 1; i++) {
        if (!units[i].next) {
            // newly added units don't have next set yet => set it
            units[i].next = units[i + 1].flowchartId;
            if (i > 0) units[i - 1].next = units[i].flowchartId;
        } else if (!flowchartIds.includes(units[i].next)) {
            // newly removed units may create broken next links => fix it
            units[i].next = units[i + 1].flowchartId;
        }
    }
    return units;
}

/**
 * @summary Apply configuration data to an object
 * @param obj {*} object / class containing methods or attributes to be set
 * @param config { functions: {}, attributes: {} } functions to call and attributes to set
 * @param callBuild {boolean} if true; call build between applying functions and attributes
 * @returns {*} updated object
 */
export function applyConfig({ obj, config = {}, callBuild = false }) {
    const { functions = {}, attributes = {} } = config;
    for (const [func, args] of Object.entries(functions)) {
        obj[func] ? args ? obj[func](args) : obj[func]() : null;
    }
    const modified = (callBuild) ? obj.build() : obj;
    for (const [key, values] of Object.entries(attributes)) {
        modified[key] = values;
    }
    return modified;
}
