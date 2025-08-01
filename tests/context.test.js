/* eslint-disable max-classes-per-file */
import ContextProvider from "@exabyte-io/ade.js/dist/js/context/ContextProvider";
import { expect } from "chai";

import { applicationContextMixin } from "../src/context/mixins/ApplicationContextMixin";
import { materialContextMixin } from "../src/context/mixins/MaterialContextMixin";
import { globalSettings } from "../src/context/providers/settings";

class MockMaterial {
    static createDefault() {
        return "defaultMockMaterial";
    }
}

class SpecificMockMaterial {
    static createDefault() {
        return "defaultSpecificMockMaterial";
    }
}

class MockApplication {
    static createDefault() {
        return "defaultMockApplication";
    }
}

class SpecificMockApplication {
    static createDefault() {
        return "defaultSpecificMockApplication";
    }
}

class ProviderEntity extends ContextProvider {
    constructor(config) {
        super(config);
        this.initApplicationContextMixin();
        this.initMaterialContextMixin();
    }
}

applicationContextMixin(ProviderEntity.prototype);
materialContextMixin(ProviderEntity.prototype);

class DerivedProviderEntity extends ProviderEntity {}

describe("Material & Application ContextMixin", () => {
    const config = { name: "test" };

    after(() => {
        globalSettings.resetDefaults();
    });

    it("uses static entity class", () => {
        globalSettings.setMaterial(MockMaterial);
        globalSettings.setApplication(MockApplication);

        const provider = new ProviderEntity(config);
        expect(provider.material).to.be.equal("defaultMockMaterial");
        expect(provider.application).to.be.equal("defaultMockApplication");
    });

    it("uses static entity class from derived class", () => {
        globalSettings.setMaterial(SpecificMockMaterial);
        globalSettings.setApplication(SpecificMockApplication);
        const provider = new DerivedProviderEntity(config);
        expect(provider.material).to.be.equal("defaultSpecificMockMaterial");
        expect(provider.application).to.be.equal("defaultSpecificMockApplication");
    });
});
