import lodash from "lodash";

import { NamedDefaultableRepetitionRuntimeItemsImportantSettingsContextAndRenderHashedInMemoryEntity } from "@exabyte-io/code.js/dist/entity";
import { getUUID } from "@exabyte-io/code.js/dist/utils";
import { UNIT_STATUSES } from "../enums";

export class BaseUnit extends NamedDefaultableRepetitionRuntimeItemsImportantSettingsContextAndRenderHashedInMemoryEntity {
    constructor(config) {
        super({
            ...config,
            status: config.status || UNIT_STATUSES.idle,
            statusTrack: config.statusTrack || [],
            flowchartId: config.flowchartId || BaseUnit.defaultFlowchartId(),
        });
    }

    static defaultFlowchartId() {
        return getUUID();
    }

    get flowchartId() {
        return this.prop("flowchartId");
    }

    get head() {
        return this.prop("head", false);
    }

    set head(bool) {
        this.setProp("head", bool);
    }

    get next() {
        return this.prop("next");
    }

    set next(flowchartId) {
        this.setProp("next", flowchartId);
    }

    get status() {
        return lodash.get(this.lastStatusUpdate, "status") || UNIT_STATUSES.idle;
    }

    set status(s) {
        this.setProp("status", s);
    }

    get lastStatusUpdate() {
        const statusTrack = this.prop("statusTrack", []).filter(s => (s.repetition || 0) === this.repetition);
        const sortedStatusTrack = lodash.sortBy(statusTrack || [], x => x.trackedAt);
        return sortedStatusTrack[sortedStatusTrack.length - 1];
    }

    get type() {
        return this.prop("type");
    }

    get isDraft() {
        return this.prop("isDraft", false);
    }

    getHashObject() {
        return { ...this.hashObjectFromRuntimeItems, type: this.type }
    }

    /**
     * Checks whether a unit is currently in a given status (e.g. idle, active, etc). The full list can be found
     * in the UNIT_STATUSES variable in enums.js.
     * @param status (String) name of the status to check
     * @returns Boolean
     */
    isInStatus(status) {
        return (this.status === status);
    }

}
