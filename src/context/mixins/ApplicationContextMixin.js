import { globalSettings } from "../providers/settings";

export function applicationContextMixin(item) {
    const properties = {
        _application: undefined,

        initApplicationContextMixin() {
            this._application =
                (this.config.context && this.config.context.application) ||
                globalSettings.Application.createDefault();
        },

        get application() {
            return this._application;
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
