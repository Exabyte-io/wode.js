export function applicationContextMixin(item) {
    const properties = {
        _application: undefined,

        initApplicationContextMixin() {
            // @ts-ignore
            if (!this.constructor.Application) {
                throw Error("ApplicationContextMixin: Application is undefined");
            }
            this._application =
                (this.config.context && this.config.context.application) ||
                this.constructor.Application.createDefault();
        },

        get application() {
            return this._application;
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
