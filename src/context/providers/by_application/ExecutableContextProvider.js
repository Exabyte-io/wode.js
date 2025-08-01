import ContextProvider from "@exabyte-io/ade.js/dist/js/context/ContextProvider";

export default class ExecutableContextProvider extends ContextProvider {
    constructor(config) {
        super({
            ...config,
            domain: "executable",
        });
    }
}
