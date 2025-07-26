import { Application } from "@exabyte-io/ade.js";
import { Made } from "@mat3ra/made";

class GlobalSettings {
    "PointsGridFormDataProvider.defaultKPPRA" = 5;

    Material = Made.Material;

    Application = Application;

    constructor() {
        this.resetDefaults();
    }

    get defaultKPPRA() {
        return this["PointsGridFormDataProvider.defaultKPPRA"];
    }

    setApplication(application) {
        this.Application = application;
    }

    setMaterial(material) {
        this.Material = material;
    }

    setDefaultKPPRA(kppra) {
        this["PointsGridFormDataProvider.defaultKPPRA"] = kppra;
    }

    resetDefaults() {
        this.Material = Made.Material;
        this.Application = Application;
        this["PointsGridFormDataProvider.defaultKPPRA"] = 5;
    }
}

export const globalSettings = new GlobalSettings();
