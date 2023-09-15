import {
    NamedDefaultableRepetitionRuntimeItemsImportantSettingsContextAndRenderHashedInMemoryEntity,
    TaggableMixin,
} from "@exabyte-io/code.js/dist/entity";
import { getUUID } from "@exabyte-io/code.js/dist/utils";
import lodash from "lodash";
import { mix } from "mixwith";

import { UNIT_STATUSES, UnitConfig } from "./types";

// eslint-disable-next-line max-len
export class BaseUnit<T extends UnitConfig> extends mix(
    NamedDefaultableRepetitionRuntimeItemsImportantSettingsContextAndRenderHashedInMemoryEntity,
).with(TaggableMixin) {
    // TODO investigate why this is needed here instead of inherited from the mixin
    public prop: <K extends keyof T>(key: K, defaultValue?: T[K]) => T[K];
    public setProp: <K extends keyof T>(key: K, value: T[K]) => void;
    private repetition: number;
    private hashObjectFromRuntimeItems: Partial<UnitConfig>;

    constructor(config: UnitConfig) {
        super({
            ...config,
            status: config.status || UNIT_STATUSES.idle,
            statusTrack: config.statusTrack || [],
            flowchartId: config.flowchartId || BaseUnit.generateFlowChartId(),
            tags: config.tags || [],
        });
    }

    static generateFlowChartId(): string {
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

    get lastStatusUpdate(): UnitConfig["statusTrack"][0] {
        const statusTrack = this.prop("statusTrack", []).filter(
            (s) => (s.repetition || 0) === this.repetition,
        );
        const sortedStatusTrack = lodash.sortBy(statusTrack || [], (x) => x.trackedAt); // TODO: check if this is the right way to sort with TS
        return sortedStatusTrack[sortedStatusTrack.length - 1];
    }

    get type() {
        return this.prop("type");
    }

    get isDraft() {
        return this.prop("isDraft", false);
    }

    getHashObject(): Partial<UnitConfig> {
        return { ...this.hashObjectFromRuntimeItems, type: this.type };
    }

    /**
     * Checks whether a unit is currently in a given status (e.g. idle, active, etc). The full list can be found
     * in the UNIT_STATUSES variable in enums.js.
     * @param status (String) name of the status to check
     * @returns Boolean
     */
    isInStatus(status) {
        return this.status === status;
    }

    clone(extraContext) {
        const flowchartIDOverrideConfigAsExtraContext = {
            flowchartId: BaseUnit.generateFlowChartId(),
            ...extraContext,
        };
        return super.clone(flowchartIDOverrideConfigAsExtraContext);
    }

}
